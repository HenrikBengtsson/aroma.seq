%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GenericDataFileSet.PARALLEL.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{dsApply.GenericDataFileSet}
\alias{dsApply.GenericDataFileSet}
\alias{GenericDataFileSet.dsApply}
\alias{dsApply,GenericDataFileSet-method}

\title{Applies a function to each file in the file set}

\description{
  Applies a function to each file in the file set.
}

\usage{
\method{dsApply}{GenericDataFileSet}(ds, FUN, ..., args=list(), skip=FALSE, verbose=FALSE, .parallel=FALSE,
  .control=list(dW = 1))
}

\arguments{
 \item{FUN}{A \code{\link[base]{function}}.}
 \item{...}{Arguments passed to \code{FUN}.}
 \item{args}{(optional) A \code{\link[base]{list}} of additional arguments
   passed to \code{FUN}.}
 \item{skip}{If \code{\link[base:logical]{TRUE}}, already processed files are skipped.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 \item{.parallel}{If \code{\link[base:logical]{TRUE}}, parallel/distributed processing is utilized.}
 \item{.control}{(internal) A named \code{\link[base]{list}} structure controlling
       the processing.}
}

\value{
  Returns nothing.
}

\seealso{
  The \pkg{BatchJobs} package is utilized for parallel/distributed
  processing.
}

\author{Henrik Bengtsson}


\keyword{internal}
\keyword{methods}
