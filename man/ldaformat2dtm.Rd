\name{ldaformat2dtm}
\alias{ldaformat2dtm}
\alias{dtm2ldaformat}
\title{Transform data from and for use with the \pkg{lda} package}
\description{
  Data from the \pkg{lda} package is transformed to a document-term
  matrix. This data format can be used to fit topic models using package
  \pkg{topicmodels}.

  Data in form of a document-term matrix is transformed to the LDA
  format used by package \pkg{lda}.
}
\usage{
ldaformat2dtm(documents, vocab)
dtm2ldaformat(x)
}
\arguments{

  \item{documents}{A \code{list} where each entry corresponds to a
    document; for each document the number of terms occurring in the
    document are stored in a \code{matrix} with two rows such that in
    each column the first entry corresponds to the vocabulary id of the
    term and the second entry to the number of times this term occurred
    in the document.}
  \item{vocab}{A \code{"character"} vector of the terms in the
    vocabulary.}
  \item{x}{An object of class \code{"DocumentTermMatrix"} as defined in
    package \pkg{tm}.}
}
\value{
  An object of class \code{"DocumentTermMatrix"} is returned by
    \code{ldaformat2dtm()} and a list with components \code{"documents"}
    and \code{"vocab"} by \code{dtm2ldaformat()}.
}
\author{
  Bettina Gruen
}
\examples{
if (require("lda")) {
  data("cora.documents", package = "lda")
  data("cora.vocab", package = "lda")
  dtm <- ldaformat2dtm(cora.documents, cora.vocab)
  cora <- dtm2ldaformat(dtm)
  all.equal(cora, list(documents = cora.documents,
                       vocab = cora.vocab))
}
}
\keyword{utilities}
