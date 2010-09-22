\name{distHellinger}
\alias{distHellinger}
\alias{distHellinger.default}
\alias{distHellinger.simple_triplet_matrix}
\title{Compute Hellinger distance}
\description{
  The Hellinger distance between the rows of two data matrices are
  determined or if the second argument is missing between the rows of
  one data matrix.
}
\usage{
\S3method{distHellinger}{default}(x, y, ...)
\S3method{distHellinger}{simple_triplet_matrix}(x, y, ...)
}
\arguments{
  \item{x}{A data matrix.}
  \item{y}{A data matrix.}
  \item{\ldots}{Currently not used.}
}
\value{
  A matrix containing the distances.
}
\author{
  Bettina Gruen
}
\keyword{cluster}
