%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  sortByName.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{sortByName}
\alias{sortByName}

\alias{sortByName.BamDataFile}
\alias{sortByName.BamDataSet}

\title{Produce sorted and indexed BAM file(s) from BAM file(s); this method will likely be superseded shortly by a SamTransform method}

\description{
  Produce sorted and indexed BAM file(s) from BAM file(s); this method will likely be superseded shortly by a SamTransform method.
}

\usage{
sortByName(...)
}

\arguments{
 \item{pathD}{The destination path.}
 \item{bIndex}{If \code{\link[base:logical]{TRUE}}, index file created after sort.}
 \item{suffix}{Filename suffix for sorted output.}
 \item{skip}{If \code{\link[base:logical]{TRUE}}, already converted files are skipped, otherwise not.}
 \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already converted files are ignored and overwritten.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link{BamDataFile}} (\code{\link{BamDataSet}}).
}

\seealso{
  Internally \code{\link[Rsamtools]{asBam}} is utilized.
}

\author{Taku Tokuyasu}


\keyword{methods}
\keyword{internal}
