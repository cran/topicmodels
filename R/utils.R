setMethod("posterior", signature(object = "TopicModel", newdata = "missing"),
function(object, newdata, ...) {
  list(terms = exp(object@beta),
       topics = object@gamma)
})

setMethod("posterior", signature(object = "TopicModel", newdata = "ANY"),
function(object, newdata, ...) {
  if (!is(newdata, "DocumentTermMatrix")) stop("newdata is of class ", dQuote(class(newdata)),
                                               " but should be of class 'DocumentTermMatrix'")
  CLASS <- strsplit(class(object), "_")[[1]]
  list(terms = exp(object@beta),
       topics = get(CLASS[1])(newdata, method = CLASS[2], 
         model = object, control = list(estimate.beta = FALSE))@gamma)
})

setGeneric("terms")
setGeneric("topics", function(x, ...) standardGeneric("topics"))

setMethod("topics", signature(x = "TopicModel"), function(x, k, threshold, ...) 
  most_likely(x, "topics", k, threshold, ...))

setMethod("terms", signature(x = "TopicModel"), function(x, k, threshold, ...) 
  most_likely(x, "terms", k, threshold, ...))

most_likely <- function(x, which = c("terms", "topics"), k, threshold, ...)
{
  which <- match.arg(which)
  labels <- if (which == "terms") x@terms else seq_len(x@k)
  post <- posterior(x)[[which]]
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
    colnames(most) <- if (which == "terms") paste("Topic", seq_len(ncol(most))) else x@documents
  } else {
    names(most) <- if (which == "terms") paste("Topic", seq_along(most)) else x@documents
  }
  return(most)
}

setGeneric("get_df", function(object, ...) standardGeneric("get_df"))

setMethod("get_df", signature(object="LDA_Gibbs"),
function(object, ...) length(object@beta))

setMethod("get_df", signature(object="LDA_VEM"),
function(object, ...) 
  as.integer(object@control@estimate.alpha) + length(object@beta))

setMethod("get_df", signature(object="CTM_VEM"),
function(object, ...) 
  (object@k - 1) * (object@k/2 + 1) + length(object@beta))

setMethod("logLik", signature(object="TopicModel"),
function(object, ...) {
  val <- sum(object@loglikelihood)
  attr(val, "df") <- get_df(object)
  attr(val, "nobs") <- object@Dim[1]
  class(val) <- "logLik"
  val
})

distHellinger <- function(x, y, ...) UseMethod("distHellinger")

distHellinger.default <- function(x, y, ...) 
{
  if (missing(y)) {
    x <- sqrt(x) 
    z <- matrix(0, nrow = nrow(x), ncol = nrow(x))
    for (k in seq_len(nrow(x)-1)) {
      z[k,-seq_len(k)] <- z[-seq_len(k),k] <- 
        sqrt(1/2 * colSums((t(x[-seq_len(k),,drop=FALSE]) - x[k,])^2))
    }
  } else {
    if(ncol(x) != ncol(y))
      stop("'x' and 'y' must have the same number of columns")
    x <- sqrt(x); y <- sqrt(y)
    z <- matrix(0, nrow = nrow(x), ncol = nrow(y))
    for (k in seq_len(nrow(y))) z[,k] <- sqrt(1/2 * colSums((t(x) - y[k,])^2))
  }
  z
}

distHellinger.simple_triplet_matrix <- function(x, y, ...) 
{
  if (!missing(y))
    stop("if x is a 'simple_triplet_matrix' y is not allowed to be specified.")
  x$v <- sqrt(x$v)
  z <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  for (k in seq_len(nrow(x)-1)) {
    x1 <- x[-seq_len(k),]
    xk <- x[k,]
    xk$i <- rep(seq_len(nrow(x1)), each = length(xk$i))
    xk$j <- rep(xk$j, nrow(x1))
    xk$v <- rep(xk$v, nrow(x1))
    xk$nrow <- nrow(x1)
    z[k,-seq_len(k)] <- z[-seq_len(k),k] <- 
      sqrt(1/2 * row_sums((xk - x1)^2))
  }
  z
}

get_terms <- function(...) terms(...)
get_topics <- function(...) topics(...)
