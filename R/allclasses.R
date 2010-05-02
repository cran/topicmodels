##**********************************************************
## control parameters

setClass("OPTcontrol",
    representation(
        iter.max = "integer",
        tol = "numeric"),
    prototype(
        iter.max = -1L,
        tol = sqrt(.Machine$double.eps)))

setClass("TopicModelcontrol",
         representation(
            verbose    = "integer",
            prefix     = "character",
            save       = "integer",
           "VIRTUAL"),
         prototype(verbose = 0L,
                   save = 0L))

setMethod("initialize", "TopicModelcontrol", function(.Object, prefix, ...) {
  if (missing(prefix)) prefix <- tempfile()
  .Object <- callNextMethod(.Object = .Object, prefix = prefix, ...)
  invisible(.Object)
})

setClass("VEMcontrol",
    representation(
        var    = "OPTcontrol",
        em     = "OPTcontrol",
        initialize = "character",
        "VIRTUAL"),
    prototype(var = new("OPTcontrol", iter.max = 500L, tol = 10^-5),
              em = new("OPTcontrol", iter.max = 1000L, tol = 10^-3),
              initialize = "random"))

setMethod("initialize", "VEMcontrol", function(.Object, initialize = "random", ...) {
  initialize <- match.arg(initialize, c("random", "seeded", "model"))
  .Object <- callNextMethod(.Object = .Object, initialize = initialize, ...)
  invisible(.Object)
})

setClass("LDAcontrol",
    representation(
        alpha = "numeric",
        "VIRTUAL"),
    contains    = "TopicModelcontrol")

setClass("LDA_VEMcontrol",
         representation(
                        estimate.alpha   = "logical"),
         contains    = c("LDAcontrol", "VEMcontrol"), 
         prototype(estimate.alpha = TRUE))

setMethod("initialize", "LDA_VEMcontrol", function(.Object, prefix, initialize = "random", ...) {
  if (missing(prefix)) prefix <- tempfile()
  .Object <- callNextMethod(.Object = .Object, prefix = prefix, ...)
  invisible(.Object)
})

setClass("LDA_Gibbscontrol",
    representation(
        delta = "numeric",
        iter = "integer"),
    contains = "LDAcontrol", 
    prototype(verbose = 0L,
              iter = 2000L))

setClass("CTM_VEMcontrol",
    representation(
                   cg = "OPTcontrol",
                   shrinkage.covariance = "logical"),
         contains = c("TopicModelcontrol", "VEMcontrol"),
         prototype(cg = new("OPTcontrol", iter.max = 500L, tol = 10^-5),
                   verbose = 5L,
                   shrinkage.covariance = FALSE))

setMethod("initialize", "CTM_VEMcontrol", function(.Object, prefix, initialize = "random", ...) {
  if (missing(prefix)) prefix <- tempfile()
  .Object <- callNextMethod(.Object = .Object, prefix = prefix, ...)
  invisible(.Object)
})
##**********************************************************
## Topic Models Objects

setClass("TopicModel",
   representation(
                  call = "call",
                  Dim = "integer",                  
                  control = "TopicModelcontrol",
                  k = "integer",
                  terms = "vector",
                  documents = "vector",
                  beta = "matrix",
                  gamma = "matrix",
                  wordassignments = "ANY",
                  "VIRTUAL"))

setClass("VEM",
         representation(
                        loglikelihood = "numeric",
                        "VIRTUAL"))

setClass("LDA",
         representation(
                        alpha = "numeric",
                        "VIRTUAL"),
         contains = "TopicModel")

setClass("LDA_VEM",
         representation(),
         contains = c("LDA", "VEM"),
         prototype(control = new("LDA_VEMcontrol")))

setClass("LDA_Gibbs",
         representation(
                        delta = "numeric"),
         contains = "LDA",
         prototype(control = new("LDA_Gibbscontrol")))

setClass("CTM",
         representation(
                        mu = "numeric",
                        Sigma = "matrix",
                        "VIRTUAL"),
         contains = "TopicModel")

setClass("CTM_VEM",
         representation(
                        nusquared = "matrix"),
         contains = c("CTM", "VEM"),
         prototype(control = new("CTM_VEMcontrol")))
         
setMethod("show", signature(object = "TopicModel"), function(object) {
  cat("A", class(object), "topic model with", object@k, "topics.\n")
})

