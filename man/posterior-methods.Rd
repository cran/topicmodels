\name{posterior-methods}
\alias{posterior,TopicModel,missing-method}
\alias{posterior,TopicModel,ANY-method}
\title{Determine posterior probabilities}
\description{
  Determine the posterior probabilities of the topics for each document and
  of the terms for each topic for a fitted topic model.} 
\usage{
\S4method{posterior}{TopicModel,missing}(object, newdata, ...)
\S4method{posterior}{TopicModel,ANY}(object, newdata, ...)
}
\arguments{
  \item{object}{An object of class "TopicModel".}
  \item{newdata}{If missing the posteriors for the original observations
    are returned.} 
  \item{\dots}{Currently not used.}
}
\author{Bettina Gruen}
\keyword{methods}


