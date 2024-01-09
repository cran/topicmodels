\name{JSS_papers}
\alias{JSS_papers}
\docType{data}
\title{JSS Papers Dublin Core Metadata}
\description{
  Dublin Core metadata for papers published in the Journal of
  Statistical Software (JSS) from 1996 until mid-2010.
}
\usage{
data("JSS_papers")
}
\format{
  A list matrix of character vectors, with rows corresponding to papers
  and the 15 columns giving the respective Dublin Core elements (variables).

  Variables \code{title} and \code{description} give the title and the
  abstract of the paper, respectively, and \code{creator} gives the
  authors (with entries character vectors with the names of the
  individual authors).
}
\details{
  Metadata were obtained from the JSS OAI repository at
  \url{https://www.jstatsoft.org/oai} via package \pkg{OAIHarvester}
  (\url{https://CRAN.R-project.org/package=OAIHarvester}).
  Records not corresponding to papers (such as book reviews) were dropped.
  
  See the documentation of package \pkg{OAIHarvester} for more
  information on Dublin Core and OAI, and
  \url{https://www.jstatsoft.org/} for information about JSS.
}
\examples{
data("JSS_papers")
## Inspect the first records:
head(JSS_papers)
## Numbers of papers by year:
table(strftime(as.Date(unlist(JSS_papers[, "date"])), "\%Y"))
## Frequent authors:
head(sort(table(unlist(JSS_papers[, "creator"])), decreasing = TRUE))
}
\keyword{datasets}
