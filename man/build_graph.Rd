\name{build_graph}
\alias{build_graph}
\title{Construct the adjacency matrix for a topic graph}
\description{
  The lasso is used to determine which edges are present in a topic
  graph.
}
\usage{
build_graph(x, lambda, and = TRUE)
}
\arguments{
  \item{x}{Object of class \code{"CTM"}.}
  \item{lambda}{Numeric in [0,1] indicates the relative bound on the
    L1-norm of the parameters.}
  \item{and}{Logical; if \code{TRUE} the graph is computed by taking the
    intersection of the neighbors, otherwise the union is determined.}
}
\value{
  Returns an adjacency matrix for the topics versus topics graph.
}
\references{
  Blei D.M., Lafferty J.D. (2007). A Correlated Topic Model of
  Science. The Annals of Applied Statistics, 1(1), 17-35. 
}
\author{David M. Blei, modified by Bettina Gruen}
\keyword{dplot}

