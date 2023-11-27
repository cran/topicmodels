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
    \item{\code{seed}:}{Object of class \code{"integer"}; used to set
      the seed in the external code for VEM estimation and to call
      \code{set.seed} for Gibbs sampling. For Gibbs sampling it can also
      be set to \code{NA} (default) to avoid changing the seed of the
      random number generator in the model fitting call.}
    \item{\code{verbose}:}{Object of class \code{"integer"}. If a
      positive integer, then the progress is reported every
      \code{verbose} iterations. If 0 (default), no output is generated
      during model fitting.}
    \item{\code{save}:}{Object of class \code{"integer"}. If a positive
      integer the estimated model is saved all \code{verbose}
      iterations. If 0 (default), no output is generated during model
      fitting.}
    \item{\code{prefix}:}{Object of class \code{"character"}; path
      indicating where to save the intermediate results.}
    \item{\code{nstart}:}{Object of class \code{"integer"}. Number of
      repeated random starts.}
    \item{\code{best}:}{Object of class \code{"logical"}; if \code{TRUE}
      only the model with the maximum (posterior) likelihood is returned,
      by default equals \code{TRUE}.}
    \item{\code{keep}:}{Object of class \code{"integer"}; if a positive
      integer, the log-likelihood is saved every \code{keep} iterations.}
    \item{\code{estimate.beta}:}{Object of class \code{"logical"};
      controls if beta, the term distribution of the topics, is fixed,
      by default equals \code{TRUE}.}
  }

  Class \code{"VEMcontrol"} contains
  \describe{
    \item{\code{var}:}{Object of class \code{"OPTcontrol"}; controls the
      variational inference for a single document, by default
      \code{iter.max} equals 500 and \code{tol} 10^-6.}
    \item{\code{em}:}{Object of class \code{"OPTcontrol"}; controls the
      variational EM algorithm, by default \code{iter.max} equals 1000
      and \code{tol} 10^-4.}
    \item{\code{initialize}:}{Object of class \code{"character"}; one of
      \code{"random"}, \code{"seeded"} and \code{"model"}, by default
      equals \code{"random"}.}
  }
  
  Class \code{"LDAcontrol"} extends class \code{"TopicModelcontrol"} and
  has the additional slots
  \describe{
    \item{\code{alpha}:}{Object of class \code{"numeric"}; initial
      value for alpha.}
  }
  
  Class \code{"LDA_VEMcontrol"} extends classes
  \code{"LDAcontrol"} and \code{"VEMcontrol"} and has the
  additional slots
  \describe{
    \item{\code{estimate.alpha}:}{Object of class \code{"logical"};
      indicates if the parameter alpha is fixed a-priori or estimated, by
      default equals \code{TRUE}.}
  }
  
  Class \code{"LDA_Gibbscontrol"} extends classes
  \code{"LDAcontrol"} and has the additional slots
  \describe{
    \item{\code{delta}:}{Object of class \code{"numeric"}; initial value
      for delta, by default equals 0.1.}
    \item{\code{iter}:}{Object of class \code{"integer"}; number of
      Gibbs iterations (after omitting the \code{burnin} iterations), by
      default equals 2000.}
    \item{\code{thin}:}{Object of class \code{"integer"}; number of
      omitted in-between Gibbs iterations, by default equals \code{iter}.}
    \item{\code{burnin}:}{Object of class \code{"integer"}; number of
      omitted Gibbs iterations at beginning, by default equals 0.}
    \item{\code{initialize}:}{Object of class \code{"character"};
      one of \code{"random"}, \code{"beta"} and \code{"z"}, by
      default equals \code{"random"}.}
  }
  
  Class \code{"CTM_VEMcontrol"} extends classes
  \code{"TopicModelcontrol"} and \code{"VEMcontrol"} and has the
  additional slots
  \describe{
    \item{\code{cg}:}{Object of class \code{"OPTcontrol"}; controls the
      conjugate gradient iterations in fitting the variational mean and
      variance per document, by default \code{iter.max} equals 500 and
      \code{tol} 10^-5.}
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
