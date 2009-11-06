\name{CTM}
\alias{CTM}
\title{Correlated Topic Model}
\description{
  Estimate a CTM model using the VEM algorithm.
}
\usage{
CTM(x, k, method = "VEM", control = NULL, model = NULL, 
    ...)
}
\arguments{
  \item{x}{Object of class \code{"DocumentTermMatrix"}}
  \item{k}{Integer; number of topics}
  \item{method}{The method to be used for fitting; currently only
    \code{method = "VEM"} is supported.}
  \item{control}{A named list of the control parameters for estimation
    or an object of class \code{"CTM_VEMcontrol"}.}
  \item{model}{Object of class \code{"CTM"} for initialization.}
  \item{\dots}{Optional arguments. Currently not used.}
}
\details{
  The C code for CTM from David M. Blei is used to estimate and fit a
  correlated topic model. 
}
\value{
  \code{CTM()} returns an object of class
  \code{"\linkS4class{CTM}"}.
}
\references{  
  Blei D.M., Lafferty J.D. (2007). A Correlated Topic Model of
  Science. The Annals of Applied Statistics, 1(1), 17-35. 
}
\author{Bettina Gruen}
\seealso{"\linkS4class{CTM_VEMcontrol}"}
\examples{
data("AssociatedPress", package = "topicmodels")
ctm <- CTM(AssociatedPress[1:20,], k = 2)
}
\keyword{models}
