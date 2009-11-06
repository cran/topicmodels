\name{TopicModelcontrol-class}
\docType{class}
\alias{OPTcontrol-class}
\alias{TopicModelcontrol-class}
\alias{CTM_VEMcontrol-class}
\alias{LDAcontrol-class}
\alias{LDA_VEMcontrol-class}
\alias{LDA_Gibbscontrol-class}
\alias{coerce,NULL,LDA_VEMcontrol-method}
\alias{coerce,list,LDA_VEMcontrol-method}
\alias{coerce,NULL,LDcontrol-method}
\alias{coerce,list,LDA_VEMcontrol-method}
\alias{coerce,NULL,CTM_VEMcontrol-method}
\alias{coerce,list,CTM_VEMcontrol-method}
\alias{coerce,NULL,OPTcontrol-method}
\alias{coerce,list,OPTcontrol-method}
\title{Different classes for controlling the estimation of topic models}
\description{
  Classes to control the estimation of topic models which are inheriting
  from the virtual base class \code{"TopicModelcontrol"}. 
}

\section{Objects from the Class}{
  Objects can be created from named lists.
}
\section{Slots}{
  Class \code{"TopicModelcontrol"} contains
  \describe{
    \item{\code{verbose}:}{Object of class \code{"integer"}; the
      estimated model is saved all \code{verbose} iterations.}
    \item{\code{prefix}:}{Object of class \code{"character"}; path
      indicating where to save the intermediate results.}
  }

  Class \code{"VEMcontrol"} contains
  \describe{
    \item{\code{var}:}{Object of class \code{"OPTcontrol"}; controls the
      variational inference for a single document.}
    \item{\code{em}:}{Object of class \code{"OPTcontrol"}; controls the
      variational EM alogorithm.}
    \item{\code{initialize}:}{Object of class \code{"character"}; one of
      \code{"random"}, \code{"seeded"} and \code{"model"}.}
  }
  
  Class \code{"LDAcontrol"} extends class \code{"TopicModelcontrol"} and has the additional
  slots
  \describe{
    \item{\code{alpha}:}{Object of class \code{"numeric"}; initial
      value for alpha.}
  }
  
  Class \code{"LDA_VEMcontrol"} extends classes
  \code{"LDAcontrol"} and \code{"VEMcontrol"} and has the
  additional slots
  \describe{
    \item{\code{estimate.alpha}:}{Object of class \code{"logical"}; by
      default equals \code{TRUE} and indicates if the parameter alpha is
      fixed a-priori or estimated.}
  }
  
  Class \code{"LDA_Gibbscontrol"} extends classes
  \code{"LDAcontrol"} and has the additional slots
  \describe{
    \item{\code{delta}:}{Object of class \code{"numeric"}; initial value
      for delta.}
    \item{\code{iter}:}{Object of class \code{"integer"}; number of
      Gibbs iterations.}
  }
  
  Class \code{"CTM_VEMcontrol"} extends classes
  \code{"TopicModelcontrol"} and \code{"VEMcontrol"} and
  has the additional slots
  \describe{
    \item{\code{cg}:}{Object of class \code{"OPTcontrol"}; controls the
      conjugate gradient iterations in fitting the variational mean and
      variance per document.}
    \item{\code{shrinkage.covariance}:}{Object of class \code{"logical"};
      indicates how the variance-covariance matrix is estimated.}
  }

  Class \code{"OPTcontrol"} contains
  \describe{
    \item{\code{iter.max}:}{Object of class \code{"integer"}; maximum
      number of iterations.}
    \item{\code{tol}:}{Object of class \code{"numeric"}; tolerance for
      convergence check.}
  }
}
\author{Bettina Gruen}
\keyword{classes}
