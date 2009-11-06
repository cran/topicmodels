\name{TopicModel-class}
\docType{class}
\alias{TopicModel-class}
\alias{LDA-class}
\alias{CTM-class}
\alias{show,TopicModel-method}
\title{Virtual class "TopicModel"}
\description{
  Fitted topic models.
}
\section{Objects from the Class}{
  Objects of class \code{"LDA"} are returned by \code{\link{LDA}()} and
  of class \code{"CTM"} by \code{\link{CTM}()}.
}
\section{Slots}{
  Class \code{"TopicModel"} contains
  \describe{
    \item{\code{call}:}{Object of class \code{"call"}}
    \item{\code{Dim}:}{Object of class \code{"integer"}; number of
      documents and terms.}
    \item{\code{control}:}{Object of class \code{"TopicModelcontrol"};
      options using for estimating the topic model.}
    \item{\code{k}:}{Object of class \code{"integer"}; number of
      topics.}
    \item{\code{terms}:}{Vector containing the terms.}
    \item{\code{documents}:}{Vector containing the document names.}
    \item{\code{beta}:}{Object of class \code{"matrix"}; parameters of
      the word distribution for each topic.}
    \item{\code{gamma}:}{Object of class \code{"matrix"}; parameters for
      the posterior topic distribution for each document.}
    \item{\code{wordassignments}:}{Object of class
      \code{"simple_triplet_matrix"}; most probable topic for each
      observed word in each document.}
  }
  
  Class \code{"VEM"} contains
  \describe{
    \item{\code{loglikelihood}:}{Object of class \code{"numeric"}; the
      log likelihood of each document given the parameters for the topic
      distribution and for the word distribution of each topic is
      approximated using the variational parameters and underestimates
      the log likelihood by the Kullback-Leibler divergence between the
      variational posterior probability and the true posterior
      probability}
  }
  
  Class \code{"LDA"} extends class \code{"TopicModel"} and has the additional
  slots
  \describe{
    \item{\code{alpha}:}{Object of class \code{"numeric"}; parameter of
      the Dirichlet distribution for topics over documents.}
  }
  
  Class \code{"LDA_Gibbs"} extends class \code{"LDA"} and has
  the additional slots
  \describe{
    \item{\code{delta}:}{Object of class \code{"numeric"}; parameter for
      the prior distribution of the words distribution for topics.}
  }
  
  Class \code{"CTM"} extends class \code{"TopicModel"} and has the additional
  slots
  \describe{
    \item{\code{mu}:}{Object of class \code{"numeric"}; mean of the
      topic distribution on the logit scale.}
    \item{\code{Sigma}:}{Object of class \code{"matrix"};
      variance-covariance matrix of topics on the logit scale.}
  }
  
  Class \code{"CTM_VEM"} extends classes \code{"CTM"} and
  \code{"VEM"} and has the additional
  slots
  \describe{
    \item{\code{nusqared}:}{Object of class \code{"matrix"}; variance of the
      variational distribution on the mu.}		
  }  
}
\author{Bettina Gruen}
\keyword{classes}
