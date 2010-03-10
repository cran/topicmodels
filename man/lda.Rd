\name{LDA}
\alias{LDA}
\title{Latent Dirichlet Allocation}
\description{
  Estimate a LDA model using the VEM algorithm or Gibbs Sampling.
}
\usage{
LDA(x, k, method = c("VEM", "Gibbs"), control = NULL,
    model = NULL, ...) 
}
\arguments{
  \item{x}{Object of class \code{"DocumentTermMatrix"}}
  \item{k}{Integer; number of topics}
  \item{method}{The method to be used for fitting; currently 
      \code{method = "VEM"} or \code{method= "Gibbs"} are
      supported.}
  \item{control}{A named list of the control parameters for estimation
    or an object of class \code{"LDAcontrol"}.}
  \item{model}{Object of class \code{"LDA"} for initialization.}
  \item{\dots}{Optional arguments. Currently not used.}
}
\details{
  The C code for LDA from David M. Blei is used to estimate and fit a
  latent dirichlet allocation model with the VEM algorithm.

  For Gibbs Sampling the C++ code from Xuan-Hieu Phan is used.
}
\value{
  \code{LDA()} returns an object of class
  \code{"\linkS4class{LDA}"}. 
}
\references{
  Blei D.M., Ng A.Y., Jordan M.I. (2003).
  Latent Dirichlet Allocation.
  \emph{Journal of Machine Learning Research}, \bold{3}, 993--1022.

  Phan X.H., Nguyen L.M., Horguchi S. (2008).
  Learning to Classify Short and Sparse Text & Web with Hidden Topics
  from Large-scale Data Collections.
  In \dQuote{Proceedings of the 17th International World Wide Web
    Conference (WWW 2008)}, pages 91--100.
  Beijing, China.
}
\author{Bettina Gruen}
\seealso{\code{"\linkS4class{LDAcontrol}"}}
\examples{
data("AssociatedPress", package = "topicmodels")
lda <- LDA(AssociatedPress[1:20,], control = list(alpha = 0.1), k = 2)
lda_inf <- LDA(AssociatedPress[21:30,], model = lda,
               control = list(em = list(iter.max = -1L)))
}
\keyword{models}
