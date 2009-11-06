CTM <- function(x, k, method = "VEM",
                control = NULL, model = NULL, ...) {
  if (!is(x, "DocumentTermMatrix")) stop("x is of class ", dQuote(class(x)))
  mycall <- match.call()
  method <- match.arg(method)
  control <- as(control, "CTM_VEMcontrol")
  if (is(model, "CTM")) control@initialize <- "model"
  if (control@initialize == "random") control@initialize <- "rand"
  else if (control@initialize == "seeded") control@initialize <- "seed" 
  if (control@initialize == "model") {
    if (!is(model, "CTM")) stop("need a model of class 'CTM' for initialization")
    x <- match_terms(x, model)
    k <- model@k
  }
  result_dir <- paste(control@prefix, "-ctm", sep = "")
  dir.create(result_dir)
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
  obj@gamma <- exp(obj@gamma)
  obj@gamma <- obj@gamma/rowSums(obj@gamma)
  new(class(obj), obj, call = mycall, control = control,
      documents = x$dimnames[[1]], terms = x$dimnames[[2]])
}
