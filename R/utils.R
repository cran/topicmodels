setMethod("posterior", signature(object = "TopicModel", newdata = "missing"),
function(object, newdata, ...) {
  list(terms = object@beta,
       topics = object@gamma)
})

setMethod("posterior", signature(object = "TopicModel", newdata = "ANY"),
function(object, newdata, ...) {
  if (!is(newdata, "DocumentTermMatrix")) stop("newdata is of class ", dQuote(class(newdata)))
  CLASS <- strsplit(class(object), "_")[[1]]
  control <- if (CLASS[2] == "VEM") list(em = list(iter.max = -1)) else list(iter = 1)
  list(terms = object@beta,
       topics = get(CLASS[1])(newdata, method = CLASS[2], 
         model = object, control = control)@gamma)
})

get_topics <- function(object, k, threshold, ...) {
  get_most_likely(object, "topics", k, threshold, ...)
}

get_terms <- function(object, k, threshold, ...) {
  get_most_likely(object, "terms", k, threshold, ...)
}

get_most_likely <- function(object, which = c("terms", "topics"), k, threshold, ...)
{
  which <- match.arg(which)
  labels <- if (which == "terms") object@terms else seq_len(object@k)
  post <- posterior(object)[[which]]
  if (missing(k)) {
    k <- ifelse(missing(threshold), 1, ncol(post))
  } else {
    k <- min(ncol(post), k)
  }
  if (!missing(threshold)) {
    most <- sapply(seq_len(nrow(post)), function(i) {
      index <- which(post[i,] > threshold)
      labels[index[index %in% order(post[i,], decreasing = TRUE)[seq_len(k)]]]
    }, ...)
  } else {
    most <- sapply(seq_len(nrow(post)), function(i) 
                    labels[order(post[i,], decreasing = TRUE)[seq_len(k)]], ...)
  }
  if (is(most, "matrix")) {
    colnames(most) <- if (which == "terms") paste("Topic", seq_len(ncol(most))) else object@documents
  } else {
    names(most) <- if (which == "terms") paste("Topic", seq_along(most)) else object@documents
  }
  return(most)
}

setGeneric("get_df", function(object, ...) standardGeneric("get_df"))

setMethod("get_df", signature(object="LDA_VEM"),
function(object, ...) {
  as.integer(object@control@estimate.alpha) + length(object@beta)
})

setMethod("get_df", signature(object="CTM_VEM"),
function(object, ...) {
  (object@k - 1) * (object@k/2 + 1) + length(object@beta)
})

setMethod("logLik", signature(object="VEM"),
function(object, ...) {
  val <- sum(object@loglikelihood)
  attr(val, "df") <- get_df(object)
  attr(val, "nobs") <- object@Dim[1]
  class(val) <- "logLik"
  val
})

