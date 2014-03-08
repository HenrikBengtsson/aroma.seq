%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  convertToSam.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{convertToSam}
\alias{convertToSam}

\alias{convertToSam.BamDataFile}
\alias{convertToSam.BamDataSet}

\title{Converts a SAM (BAM) file (set of files) set into a BAM (SAM) file (set of files)}

\description{
  Converts a SAM (BAM) file (set of files) set into a BAM (SAM) file (set of files).
  When converting to BAM, the generated BAM files are all sorted and indexed as well.
}

\usage{
convertToSam(...)
}

\arguments{
 \item{path}{The destination path.}
 \item{skip}{If \code{\link[base:logical]{TRUE}}, already converted files are skipped, otherwise not.}
 \item{overwrite}{If \code{\link[base:logical]{TRUE}}, already converted files are ignored and overwritten.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 \item{...}{Not used.}
}

\value{
  \code{convertToBam()} returns a (sorted and indexed) \code{\link{BamDataFile}} (\code{\link{BamDataSet}}).
  \code{convertToSam()} returns a \code{\link{SamDataFile}} (\code{\link{SamDataSet}}).
}

\seealso{
  Internally \code{\link[Rsamtools]{asBam}} or \code{\link[Rsamtools]{asSam}} is utilized.
}

\author{Henrik Bengtsson}


\keyword{methods}
\keyword{internal}
