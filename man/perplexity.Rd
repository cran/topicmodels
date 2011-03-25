\name{perplexity}
\alias{perplexity}
\alias{perplexity,VEM,missing-method}
\alias{perplexity,ANY,matrix-method}
\alias{perplexity,ANY,DocumentTermMatrix-method}
\alias{perplexity,VEM,simple_triplet_matrix-method}
\alias{perplexity,Gibbs,simple_triplet_matrix-method}
\alias{perplexity,Gibbs_list,simple_triplet_matrix-method}
\alias{perplexity,list,missing-method}
\alias{perplexity,list,simple_triplet_matrix-method}
\title{Methods for Function perplexity}
\description{
  Determine the perplexity of a fitted model.
}
\usage{
perplexity(object, newdata, \ldots)

\S4method{perplexity}{VEM,simple_triplet_matrix}(object, newdata, control, \ldots)

\S4method{perplexity}{Gibbs,simple_triplet_matrix}(object, newdata, control, use_theta = TRUE,
estimate_theta = TRUE, \ldots)

\S4method{perplexity}{Gibbs_list,simple_triplet_matrix}(object, newdata, control, use_theta  = TRUE,
estimate_theta = TRUE, \ldots)

}
\arguments{
  \item{object}{Object of class \code{"TopicModel"} or \code{"Gibbs_list"}.}
  \item{newdata}{If missing, the perplexity for the data to which the
    model was fitted is determined. For objects fitted using Gibbs sampling
    \code{newdata} needs to be specified.}
  \item{control}{If missing, the \code{control} of the fitted model is
    used with suitable changes of the relevant parameters (see
    \bold{Details}).}
  \item{use_theta}{Object of class \code{"logical"}. If \code{TRUE}
    the estimated topic distributions for the documents are
    used. Otherwise equal weights are assigned to the topics for each document.}
  \item{estimate_theta}{Object of class \code{"logical"}. If \code{FALSE} the
    data provided is assumed to be the same as the data used for fitting the
    model. The topic distributions therefore do not need to be estimated
    and the data in \code{newdata} is used for weighting the
    term-document occurrences.}
  \item{\dots}{Further arguments passed to the different methods.}
}
\details{
  The specified control is modified to ensure that (1)
  \code{estimate.beta=FALSE} and (2) \code{nstart=1}.

  For \code{"Gibbs_list"} objects the \code{control} is further modified
  to have (1) \code{iter=thin} and (2) \code{best=TRUE} and the model is
  fitted to the new data with this control for each available
  iteration. The perplexity is then determined by averaging over the
  same number of iterations.

  If a \code{list} is supplied as \code{object}, it is assumed that it
  consists of several models which were fitted using different starting
  configurations.
}
\value{
  A numeric value.
}
\references{
  Blei D.M., Ng A.Y., Jordan M.I. (2003).
  Latent Dirichlet Allocation.
  \emph{Journal of Machine Learning Research}, \bold{3}, 993--1022.

  Griffiths T.L., Steyvers, M. (2004).
  Finding Scientific Topics.
  \emph{Proceedings of the National Academy of Sciences of the
    United States of America}, \bold{101}, Suppl. 1, 5228--5235.
  
  Newman D., Asuncion A., Smyth P., Welling M. (2009).
  Distributed Algorithms for Topic Models.
  \emph{Journal of Machine Learning Research}, \bold{10}, 1801--1828.
}
\author{Bettina Gruen}
\keyword{methods}
