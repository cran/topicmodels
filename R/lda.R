match_terms <- function(x, model) {
  if (!(is.null(model@terms) | is.null(colnames(x)))) {
    x <- x[,which(colnames(x) %in% model@terms)]
    x <- x[,order(match(colnames(x), model@terms))]
    js <- match(colnames(x), model@terms)
    x$j <- js[x$j]
    x$ncol <- model@Dim[2]
    dimnames(x)[[2]] <- model@terms
  } else if (ncol(x) != model@Dim[2]) stop("the number of terms in the input matrix and the fitted model need to match")
  x
}

LDA_registry <- list(LDA_VEM.fit = c("VEM", "LDA_VEM", "LDA_VEM.fit"),
                 LDA_Gibbs.fit = c("Gibbs", "LDA_Gibbs", "LDA_Gibbs.fit"))

LDA <- function(x, k, method = "VEM", control = NULL, model = NULL, ...)
{
  if (is(x, "DocumentTermMatrix")) {
    if (!any(attr(x, "Weighting") %in% c("term frequency", "tf"))) {
      stop("\nDocumentTermMatrix needs to have a term frequency weighting")
    }
  } else if (!is(x, "simple_triplet_matrix")) {
    x <- slam::as.simple_triplet_matrix(x)
  }
  if (!all.equal(x$v, as.integer(x$v)))
    stop("\nInput matrix needs to contain integer entries")
  if (!all(slam::row_sums(x) > 0)) stop("\nEach row of the input matrix needs to contain at least one non-zero entry")
  mycall <- match.call()

  if (!is.null(model)) {
    x <- match_terms(x, model)
    k <- model@k
  }

  if (as.integer(k) != k || as.integer(k) < 2) stop("\nk needs to be an integer of at least 2")

  if(missing(method) && !missing(model))
    method <- paste(class(model), "fit", sep = ".")
  if(!is.function(method)) {
    MATCH <- which(sapply(LDA_registry, function(x) length(grep(tolower(method), tolower(x)))) > 0)
    if (!length(MATCH) == 1)
      stop("\nMethod not specified correctly")
    method <- get(names(LDA_registry)[MATCH])
  }

  method(x, k, control, model, mycall, ...)
}

LDA_VEM.fit <- function(x, k, control = NULL, model = NULL, call, ...) {
  control <- as(control, "LDA_VEMcontrol")
  if (length(control@seed) != control@nstart)
    stop(paste("\nneed ", control@nstart, " seeds", sep = ""))
  if (length(control@alpha) == 0)  {
    control@alpha <- if (!is.null(model)) model@alpha else 50/k
  }
  if (is.null(model)) {
    if (control@initialize == "model")
      stop("\nNeed a model of class 'LDA_VEM' for initialization")
  }
  else control@initialize <- "model"
  if (!control@estimate.beta) control@em@iter.max <- -1L
  result_dir <- path.expand(paste(control@prefix, "-lda", sep = ""))

  if (control@save) dir.create(result_dir, showWarnings = FALSE)

  obj <- vector("list", control@nstart)
  for (i in seq_len(control@nstart)) {
    control_i <- control
    control_i@seed <- control@seed[i]
    obj[[i]] <- .Call("rlda", 
                      ## simple_triplet_matrix
                      as.integer(x$i),
                      as.integer(x$j),
                      as.numeric(x$v),
                      as.integer(x$nrow),
                      as.integer(x$ncol),                 
                      ## LDAcontrol
                      control_i,
                      ## number of topics
                      as.integer(k),
                      ## directory for output files
                      result_dir, 
                      ## initial model
                      model,
                      PACKAGE = "topicmodels")
    obj[[i]]@gamma <- obj[[i]]@gamma/rowSums(obj[[i]]@gamma)
    obj[[i]] <- new(class(obj[[i]]), obj[[i]], call = call, control = control_i,
                    documents = x$dimnames[[1]], terms = x$dimnames[[2]], n = as.integer(sum(x$v)))
  }
  if (control@best) {
    MAX <- which.max(sapply(obj, logLik))
    if (length(MAX)) {
      obj <- obj[[MAX]]
    } else warning("problem selecting best fitting model")
  }
  obj
}

LDA_Gibbs.fit <- function(x, k, control = NULL, model = NULL, call, ...) {
  if (!is.null(model) && is(control, "list") && !"delta" %in% names(control)) control <- c(control, delta = model@delta)
  control <- as(control, "LDA_Gibbscontrol")
  if (length(control@seed) != control@nstart)
    stop(paste("\nneed ", control@nstart, " seeds", sep = ""))
  if (length(control@alpha) == 0)  {
    control@alpha <- if (!is.null(model)) model@alpha else 50/k
  }
  
  result_dir <- path.expand(paste(control@prefix, "-lda", sep = ""))
  if (control@save) dir.create(result_dir, showWarnings = FALSE)
  CONTROL_i <- control
  CONTROL_i@iter <- control@burnin + control@thin
  obj <- vector("list", control@nstart)
  for (i in seq_len(control@nstart)) {
    CONTROL_i@seed <- CONTROL_i@seed[i]
    obj[[i]] <- list(.Call("rGibbslda", 
                           ## simple_triplet_matrix
                           as.integer(x$i),
                           as.integer(x$j),
                           as.numeric(x$v),
                           as.integer(x$nrow),
                           as.integer(x$ncol),                 
                           ## LDAcontrol
                           CONTROL_i,
                           ## initialize
                           is.null(model),
                           ## number of topics
                           as.integer(k),
                           ## directory for output files
                           result_dir,
                           ## initial model
                           model,
                           PACKAGE = "topicmodels"))
    obj[[i]][[1]] <- new(class(obj[[i]][[1]]), obj[[i]][[1]], call = call, control = CONTROL_i,
                         documents = x$dimnames[[1]], terms = x$dimnames[[2]], n = as.integer(sum(x$v)))
    iterations <- unique(c(seq(CONTROL_i@iter, control@burnin + control@iter, by = control@thin),
                           control@burnin + control@iter))
    if (length(iterations) > 1) {
      for (j in seq_along(iterations)[-1]) {
        CONTROL_i@iter <- diff(iterations)[j-1]
        obj[[i]][[j]] <- .Call("rGibbslda", 
                               ## simple_triplet_matrix
                               as.integer(x$i),
                               as.integer(x$j),
                               as.numeric(x$v),
                               as.integer(x$nrow),
                               as.integer(x$ncol),                 
                               ## LDAcontrol
                               CONTROL_i,
                               ## initialize
                               FALSE,
                               ## number of topics
                               as.integer(k),
                               ## directory for output files
                               result_dir, 
                               ## initial model
                               obj[[i]][[j-1]],
                               PACKAGE = "topicmodels")
        obj[[i]][[j]] <- new(class(obj[[i]][[j]]), obj[[i]][[j]], call = call, control = CONTROL_i,
                             documents = x$dimnames[[1]], terms = x$dimnames[[2]], n = as.integer(sum(x$v)))
      }
      if (control@best) obj[[i]] <- obj[[i]][[which.max(sapply(obj[[i]], logLik))]]
    } else if (control@best) obj[[i]] <- obj[[i]][[1]]
  }    
  if (control@best) {
    MAX <- which.max(sapply(obj, logLik))
    if (length(MAX)) {
      obj <- obj[[MAX]]
    } else warning("no finite likelihood")
  } else {
    obj <- lapply(obj, function(x) new("Gibbs_list", fitted = x))
    if (control@nstart == 1) obj <- obj[[1]]
  }
  obj
}
