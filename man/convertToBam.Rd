%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  convertToBam.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{convertToBam}
\alias{convertToBam}

\alias{convertToBam.SamDataFile}
\alias{convertToBam.SamDataSet}

\title{Converts a SAM file (set of files) set into a (sorted and indexed) BAM file (set of files)}

\description{
  Converts a SAM file (set of files) set into a (sorted and indexed) BAM file (set of files).
}

\usage{
convertToBam(...)
}

\arguments{
 \item{path}{The destination path.}
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

\author{Henrik Bengtsson}


\keyword{methods}
\keyword{internal}
