\name{build_graph}
\alias{build_graph}
\title{Construct the adjacency matrix for a topic graph}
\description{ 
  The lasso is used to determine which edges are present in a topic
  graph. The original R code was written by David M. Blei and co-authors
  and is available together with the C code for fitting the CTM.
}
\usage{
build_graph(x, lambda, and = TRUE)
}
\arguments{
  \item{x}{Object of class \code{"CTM"}.}
  \item{lambda}{Numeric in \eqn{[0,1]} indicating the relative bound on
    the L1-norm of the parameters.}
  \item{and}{Logical; if \code{TRUE} the graph is computed by taking the
    intersection of the neighbors, otherwise the union is determined.}
}
\value{
  Returns an adjacency matrix for the topics versus topics graph.
}
\references{
  Blei D.M., Lafferty J.D. (2007).
  A Correlated Topic Model of Science.
  \emph{The Annals of Applied Statistics}, \bold{1}(1), 17--35. 
}
\author{David M. Blei, modified by Bettina Gruen.}
\keyword{dplot}
