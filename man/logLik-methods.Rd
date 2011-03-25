\name{logLik-methods}
\docType{methods}
\title{Methods for Function logLik}
\alias{logLik,TopicModel-method}
\alias{logLik,Gibbs_list-method}
\description{Compute the log-likelihood.}
\section{Methods}{
  \describe{
    \item{object = TopicModel:}{Compute the log-likelihood of a
      \code{"TopicModel"} object. For \code{"VEM"} objects the sum of
      the log-likelihood of all documents given the parameters for the
      topic distribution and for the word distribution of each topic is
      approximated using the variational parameters and underestimates
      the log-likelihood by the Kullback-Leibler divergence between the
      variational posterior probability and the true posterior
      probability.}
    \item{object = Gibbs_list:}{Compute the log-likelihoods of the
      \code{"TopicModel"} objects contained in the \code{"Gibbs_list"}
      object.}
  }
}
\keyword{methods}
