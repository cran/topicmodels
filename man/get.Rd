\name{get}
\Rdversion{1.1}
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
get_terms(object, k, threshold, ...)
get_topics(object, k, threshold, ...)
}
\arguments{
  \item{object}{Object of class \code{"TopicModel"}.}
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
