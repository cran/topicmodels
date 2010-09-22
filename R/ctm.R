CTM_registry <- list(CTM_VEM.fit = c("VEM", "CTM_VEM", "CTM_VEM.fit"))

CTM <- function(x, k, method = "VEM", control = NULL, model = NULL, ...) {
  if (!is(x, "DocumentTermMatrix")) stop("\nx is of class ", dQuote(class(x)))
  if (!all(slam::row_sums(x) > 0)) stop("\nAll documents in the DocumentTermMatrix need to contain at least one term")
  if (!all.equal(x$v, as.integer(x$v)))
    stop("\nDocumentTermMatrix needs to contain integer entries")
  mycall <- match.call()
  
  if (!is.null(model)) {
    x <- match_terms(x, model)
    k <- model@k
  }

  if (as.integer(k) != k || as.integer(k) < 2) stop("k needs to be an integer of at least 2")  

  if(is.null(method))
    method <- if (!missing(model)) paste(class(model), "fit", sep = ".") else CTM_VEM.fit
  if(!is.function(method)) {
    MATCH <- which(sapply(CTM_registry, function(x) length(grep(tolower(method), tolower(x)))) > 0)
    if (!length(MATCH) == 1)
      stop("\nMethod not specified correctly")
    method <- get(names(CTM_registry)[MATCH])
  }

  method(x, k, control, model, mycall, ...)
}

CTM_VEM.fit <- function(x, k, control = NULL, model = NULL, call, ...) {
  control <- as(control, "CTM_VEMcontrol")
  if (control@initialize == "random") control@initialize <- "rand"
  else if (control@initialize == "seeded") control@initialize <- "seed" 
  else if (control@initialize == "model") {
    if (!is(model, "CTM")) stop("need a model of class 'CTM' for initialization")
  }
  if (is(model, "CTM")) control@initialize <- "model"
  result_dir <- path.expand(paste(control@prefix, "-ctm", sep = ""))
  if (control@save) dir.create(result_dir, showWarnings = FALSE)
  if (!control@estimate.beta) control@em@iter.max <- -1L
  obj <- .Call("rctm",
               ## simple_triplet_matrix
               as.integer(x$i),
               as.integer(x$j),
               as.numeric(x$v),
               as.integer(x$nrow),
               as.integer(x$ncol),                 
               ## CTMcontrol
               control,
               ## number of topics
               as.integer(k),
               ## directory for output files
               result_dir,
               ## initial model
               model,
               PACKAGE = "topicmodels")
  obj@gamma <- cbind(exp(obj@gamma), 1)
  obj@gamma <- obj@gamma/rowSums(obj@gamma)
  new(class(obj), obj, call = call, control = control,
      documents = x$dimnames[[1]], terms = x$dimnames[[2]])
}
