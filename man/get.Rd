\name{terms/topics}
\docType{methods}
\alias{topics}
\alias{topics,TopicModel-method}
\alias{terms,TopicModel-method}
\alias{get_terms}
\alias{get_topics}
\title{
  Extract most likely terms or topics.
}
\description{
  Function to extract the most likely terms for each topic or the most
  likely topics for each document.
}
\usage{
\S4method{terms}{TopicModel}(x, k, threshold, ...)
\S4method{topics}{TopicModel}(x, k, threshold, ...)
}
\arguments{
  \item{x}{Object of class \code{"TopicModel"}.}
  \item{k}{The maximum number of terms/topics returned. By default
    set to 1 if no threshold is given.}
  \item{threshold}{Only the terms/topics which are more likely than
    the threshold are returned.}
  \item{\dots}{Further arguments passed to \code{sapply}.}
}
\value{
  A list or matrix containing the most likely terms for each topic or
  the most likely topics for each document.
}
\author{
  Bettina Gruen
}
\seealso{\code{\link[topicmodels]{posterior-methods}}}
\keyword{utilities}
