\name{NEWS}
\title{News for Package 'topicmodels'}
\encoding{UTF-8}
\section{Changes in topicmodels version 0.1-9}{
  \itemize{
    \item A bug in function \code{posterior} was fixed where the
    rownames of the wrong object were used. Thanks to Benjamin S. Porter
    for pointing the problem out.
    \item Dependency structure changed such that some packages are now
    only imported.
    \item The information printed during the VEM algorithm when
    \code{verbose} is larger than 0 was improved.
  }
}
\section{Changes in topicmodels version 0.1-8}{
  \itemize{
    \item The code in the vignette for removing HTML markup was modified
    due to changes in package \pkg{XML}.
  }
}
\section{Changes in topicmodels version 0.1-7}{
  \itemize{
    \item A memory leak in the code of the fit function for LDA with method
      \code{"VEM"} was corrected. Thanks to Ramis Yamilov for pointing
      the problem out. 
  }
}
\section{Changes in topicmodels version 0.1-6}{
  \itemize{
    \item The included dataset AssociatedPress had row names which were of
      type integer and not of type character. The object was re-saved
      omitting the row names.
  }
}
\section{Changes in topicmodels version 0.1-5}{
  \itemize{
    \item Vignettes moved from /inst/doc to /vignettes.
    \item The source code for fitting the model using Gibbs sampling was
      modified because the code did not compile on Solaris. Thanks to
      Prof. Brian D. Ripley for pointing the problem out.
    \item \code{dtm2ldaformat()} was modified to ensure that the resulting matrices
      for the documents contain integers. In addition \code{dtm2ldaformat()}
      and \code{ldaformat2dtm()} were changed to also work for document-term
      matrices containing empty documents and an argument was introduced
      to indicate if empty documents should be removed. Thanks to Eu Jin
      Lok for pointing the problems out.
  }
}
\section{Changes in topicmodels version 0.1-4}{
  \itemize{
    \item Missing 'Suggests' entries added in the DESCRIPTION file. Thanks
      to Prof. Brian D. Ripley for pointing the problem out.
  }
}
\section{Changes in topicmodels version 0.1-3}{
  \itemize{
    \item Name tags for Rd files changed to not contain slashes. Thanks to
      Prof. Brian D. Ripley for pointing the problem out as indicated in
      bug PR14707.
  }
}
\section{Changes in topicmodels version 0.1-2}{
  \itemize{
    \item A small bug fixed when saving interim results for fitting a LDA
      model using Gibbs sampling. Thanks to Nicholas Switanek for
      pointing the problem out.
  }
}
\section{Changes in topicmodels version 0.1-1}{
  \itemize{
    \item Makevars.win changed due to changes on CRAN for making libgsl for
      Windows. Thanks to Prof. Brian D. Ripley for pointing that out.
  }
}
\section{Changes in topicmodels version 0.1-0}{
  \itemize{
    \item The package vignette has been published in the Journal of
      Statistical Software, Volume 40, Issue 13
      (\url{http://www.jstatsoft.org/v40/i13}), and the paper should be
      used as citation for the package, run
      \code{citation("topicmodels")} for details. 
  }
}
\section{Changes in topicmodels version 0.0-11}{
  \itemize{
    \item C code changed to allow the package to compile on Solaris
      systems. Thanks to Prof. Brian D. Ripley for pointing the problems
      out and recommending suitable changes.
  }
}
\section{Changes in topicmodels version 0.0-10}{
  \itemize{
    \item C code changed to avoid warnings of unused variables.
  }
}
\section{Changes in topicmodels version 0.0-9}{
  \itemize{
    \item The slots for documents and terms names are not restricted to be
      of class \code{"vector"} any more to allow for document-term matrices
      where no row and/or column names are provided.
  }
}
\section{Changes in topicmodels version 0.0-8}{
  \itemize{
    \item A function \code{perplexity()} added for model validation and selection. 
    \item The input data for \code{LDA()} and \code{CTM()} can now either be a
      \code{"DocumentTermMatrix"} with term-frequency weighting or an object
      coercible to a \code{"simple_triplet_matrix"} with integer entries.
    \item A bug in the C++ Gibbs sampling code fixed for the random number
      generation. Thanks to Uwe Ligges for pointing the problem out
      which he noted when checking the package for the Windows
      platform.
    \item New control arguments added for keeping intermediate
      log-likelihood values during estimation and running repeated runs
      with random initilization.  In addition the number of iterations
      made is now saved with the fitted model.
    \item Functions \code{ldaformat2dtm()} and \code{dtm2ldaformat()} added to transform data
    from the \pkg{lda} package into a \code{"DocumentTermMatrix"} object and vice
      versa.
    \item Bug fixed in rctm.c where for \code{estimate.beta = FALSE} one EM step
      was performed.
  }
}
\section{Changes in topicmodels version 0.0-7}{
  \itemize{
    \item The control for topic models now also has a \code{seed}
      argument to ensure reproducibility of results and a
      \code{estimate.beta} argument which can be used to fix the term
      distribution over topics after initialization.
    \item The control for Gibbs sampling allows to specify to return
      repeated draws in a list using arguments \code{burnin},
      \code{thin} and \code{iter}. 
    \item In slot beta for class \code{"TopicModel"} the log parameters are stored
      to have a higher accuracy for the VEM code if parameter values are
      close to zero.
    \item Call to assert removed in C code to avoid termination of \R.
    \item Class \code{"TopicModel"} now has a slot
      \code{loglikelihood}. For models fitted using Gibbs sampling this
      contains the loglikelihood of the corpus, for VEM fitted models
      the vector of loglikelihoods for each document separately.
  }
}
\section{Changes in topicmodels version 0.0-6}{
  \itemize{
    \item Memory bug fixed in \code{returnObjectGibbsLDA}.
    \item A slot \code{save} is added to the control objects to specify if the
      results and with which step size intermediate results are saved
      into files.
  }
}
\section{Changes in topicmodels version 0.0-5}{
  \itemize{
    \item Header files changed in utilities.cpp following an advice by
      Prof. Brian D. Ripley.
  }
}
\section{Changes in topicmodels version 0.0-4}{
  \itemize{
    \item Code for installing the package \pkg{corpus.JSS.papers} in the
      vignette improved.
    \item \code{dir.create()} now called with \code{showWarnings = FALSE}.
    \item Bug fixed in \code{get_most_likely()} for maximum possible k.
    \item First version released on CRAN: 0.0-3.
  }
}
