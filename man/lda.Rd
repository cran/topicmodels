\name{LDA}
\alias{LDA}
\title{Latent Dirichlet Allocation}
\description{
  Estimate a LDA model using for example the VEM algorithm or Gibbs Sampling.
}
\usage{
LDA(x, k, method = "VEM", control = NULL, model = NULL, ...) 
}
\arguments{
  \item{x}{Object of class \code{"DocumentTermMatrix"} with
    term-frequency weighting or an object coercible to a
    \code{"simple_triplet_matrix"} with integer entries.}
  \item{k}{Integer; number of topics.}
  \item{method}{The method to be used for fitting; currently 
      \code{method = "VEM"} or \code{method= "Gibbs"} are
      supported.}
  \item{control}{A named list of the control parameters for estimation
    or an object of class \code{"LDAcontrol"}.}
  \item{model}{Object of class \code{"LDA"} for initialization.}
  \item{\dots}{Optional arguments. For \code{method = "Gibbs"} an
    additional argument \code{seedwords} can be specified as a matrix or
    an object of class \code{"simple_triplet_matrix"}; the default is
    \code{NULL}.}
}
\details{
  The C code for LDA from David M. Blei and co-authors is used to
  estimate and fit a latent dirichlet allocation model with the VEM
  algorithm. For Gibbs Sampling the C++ code from Xuan-Hieu Phan and
  co-authors is used.

  When Gibbs sampling is used for fitting the model, seed words with
  their additional weights for the prior parameters can be specified in
  order to be able to fit seeded topic models.
}
\value{
  \code{LDA()} returns an object of class \code{"\linkS4class{LDA}"}. 
}
\references{
  Blei D.M., Ng A.Y., Jordan M.I. (2003).
  Latent Dirichlet Allocation.
  \emph{Journal of Machine Learning Research}, \bold{3}, 993--1022.

  Phan X.H., Nguyen L.M., Horguchi S. (2008).
  Learning to Classify Short and Sparse Text & Web with Hidden Topics
  from Large-scale Data Collections.
  In Proceedings of the 17th International World Wide Web Conference
  (WWW 2008), pages 91--100, Beijing, China.

  Lu, B., Ott, M., Cardie, C., Tsou, B.K. (2011).
  Multi-aspect Sentiment Analysis with Topic Models.
  In Proceedings of the 2011 IEEE 11th International Conference on Data
  Mining Workshops, pages 81--88.
}
\author{Bettina Gruen}
\seealso{\code{"\linkS4class{LDAcontrol}"}}
\examples{
data("AssociatedPress", package = "topicmodels")
lda <- LDA(AssociatedPress[1:20,], control = list(alpha = 0.1), k = 2)
lda_inf <- posterior(lda, AssociatedPress[21:30,])
}
\keyword{models}
