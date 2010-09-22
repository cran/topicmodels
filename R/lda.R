match_terms <- function(x, model) {
  x <- x[,which(colnames(x) %in% model@terms)]
  x <- x[,order(match(colnames(x), model@terms))]
  js <- match(colnames(x), model@terms)
  x$j <- js[x$j]
  x$ncol <- model@Dim[2]
  dimnames(x)[[2]] <- model@terms
  x
}

LDA_registry <- list(LDA_VEM.fit = c("VEM", "LDA_VEM", "LDA_VEM.fit"),
                 LDA_Gibbs.fit = c("Gibbs", "LDA_Gibbs", "LDA_Gibbs.fit"))

LDA <- function(x, k, method = "VEM", control = NULL, model = NULL, ...)
{
  if (!is(x, "DocumentTermMatrix")) stop("\nx is of class ", dQuote(class(x)))
  if (!all(slam::row_sums(x) > 0)) stop("\nAll documents in the DocumentTermMatrix need to contain at least one term")
  if (!all.equal(x$v, as.integer(x$v)))
    stop("\nDocumentTermMatrix needs to contain integer entries")
  mycall <- match.call()
  
  if (!is.null(model)) {
    x <- match_terms(x, model)
    k <- model@k
  }

  if (as.integer(k) != k || as.integer(k) < 2) stop("\nk needs to be an integer of at least 2")

  if(is.null(method))
    method <- if (!missing(model)) paste(class(model), "fit", sep = ".") else LDA_VEM.fit
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
  if (length(control@alpha) == 0)  {
    control@alpha <- if (!is.null(model)) model@alpha else 50/k
  }
  if (is.null(model)) {
    if (control@initialize == "model")
      stop(paste("\nNeed a model of class 'LDA_VEM' for initialization", sep = ""))
  }
  else control@initialize <- "model"
  if (!control@estimate.beta) control@em@iter.max <- -1L
  result_dir <- path.expand(paste(control@prefix, "-lda", sep = ""))
  if (control@save) dir.create(result_dir, showWarnings = FALSE)
  obj <- .Call("rlda", 
               ## simple_triplet_matrix
               as.integer(x$i),
               as.integer(x$j),
               as.numeric(x$v),
               as.integer(x$nrow),
               as.integer(x$ncol),                 
               ## LDAcontrol
               control,
               ## number of topics
               as.integer(k),
               ## directory for output files
               result_dir, 
               ## initial model
               model,
               PACKAGE = "topicmodels")
  obj@gamma <- obj@gamma/rowSums(obj@gamma)
  new(class(obj), obj, call = call, control = control,
      documents = x$dimnames[[1]], terms = x$dimnames[[2]])
}

LDA_Gibbs.fit <- function(x, k, control = NULL, model = NULL, call, ...) {
  if (!is.null(model) && !"delta" %in% names(control)) control <- c(control, delta = model@delta)
  control <- as(control, "LDA_Gibbscontrol")
  if (length(control@alpha) == 0)  {
    control@alpha <- if (!is.null(model)) model@alpha else 50/k
  }
  
  result_dir <- path.expand(paste(control@prefix, "-lda", sep = ""))
  if (control@save) dir.create(result_dir, showWarnings = FALSE)
  CONTROL <- control
  CONTROL@iter <- control@burnin + control@thin
  obj <- .Call("rGibbslda", 
               ## simple_triplet_matrix
               as.integer(x$i),
               as.integer(x$j),
               as.numeric(x$v),
               as.integer(x$nrow),
               as.integer(x$ncol),                 
               ## LDAcontrol
               CONTROL,
               ## initialize
               is.null(model),
               ## number of topics
               as.integer(k),
               ## directory for output files
               dir,
               ## initial model
               model,
               PACKAGE = "topicmodels")
  obj <- new(class(obj), obj, call = call, control = CONTROL,
             documents = x$dimnames[[1]], terms = x$dimnames[[2]])
  iterations <- unique(c(seq(CONTROL@iter, control@burnin + control@iter, by = control@thin),
                         control@burnin + control@iter))
  if (length(iterations) > 1) {
    obj <- list(obj)
    for (i in seq_along(iterations)[-1]) {
      CONTROL@iter <- diff(iterations)[i-1]
      obj[[i]] <- .Call("rGibbslda", 
                        ## simple_triplet_matrix
                        as.integer(x$i),
                        as.integer(x$j),
                        as.numeric(x$v),
                        as.integer(x$nrow),
                        as.integer(x$ncol),                 
                        ## LDAcontrol
                        CONTROL,
                        ## initialize
                        FALSE,
                        ## number of topics
                        as.integer(k),
                        ## directory for output files
                        result_dir, 
                        ## initial model
                        obj[[i-1]],
                        PACKAGE = "topicmodels")
      obj[[i]] <- new(class(obj[[i]]), obj[[i]], call = call, control = CONTROL,
                      documents = x$dimnames[[1]], terms = x$dimnames[[2]])
    }
    if (control@best) obj <- obj[[which.max(sapply(obj, logLik))]]
  }
  obj
}
