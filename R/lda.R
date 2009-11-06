match_terms <- function(x, model) {
  x <- x[,which(colnames(x) %in% model@terms)]
  x <- x[,order(match(colnames(x), model@terms))]
  js <- match(colnames(x), model@terms)
  x$j <- js[x$j]
  x$ncol <- model@Dim[2]
  dimnames(x)[[2]] <- model@terms
  x
}

LDA <- function(x, k, method = c("VEM", "Gibbs"), 
                control = NULL, model = NULL, ...) {
  if (!is(x, "DocumentTermMatrix")) stop("x is of class ", dQuote(class(x)))
  mycall <- match.call()
  method <- match.arg(method)
  control <- as(control, paste("LDA_", method, "control", sep = ""))
  if (!is.null(model)) {
    if (!is(model, paste("LDA_", method, sep = ""))) stop(paste("need a model of class 'LDA_", method, "' for initialization", sep = ""))
    x <- match_terms(x, model)
    control@alpha <- model@alpha
    if (method == "Gibbs") control@delta <- model@delta
    k <- model@k
  }
  if (method == "VEM") {
    if (is(model, paste("LDA_", method, sep = ""))) control@initialize <- "model"
    if (control@initialize == "model") {
      if (is.null(model)) stop(paste("need a model of class 'LDA_", method, "' for initialization", sep = ""))
    }
  }
  result_dir <- paste(control@prefix, "-lda", sep = "")
  dir.create(result_dir)
  if (length(control@alpha) == 0) control@alpha <- 50/k
  if (method == "VEM") {
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
  } else {
    if (length(control@delta) == 0) control@delta <- 0.1
    obj <- .Call("rGibbslda", 
                 ## simple_triplet_matrix
                 as.integer(x$i),
                 as.integer(x$j),
                 as.numeric(x$v),
                 as.integer(x$nrow),
                 as.integer(x$ncol),                 
                 ## LDAcontrol
                 control,
                 ## character: seeded, random, model
                 initialize,
                 ## number of topics
                 as.integer(k),
                 ## directory for output files
                 result_dir,
                 ## initial model
                 model,
                 PACKAGE = "topicmodels")
  }
  new(class(obj), obj, call = mycall, control = control,
      documents = x$dimnames[[1]], terms = x$dimnames[[2]])
}
