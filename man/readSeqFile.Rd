%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  readSeqFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readSeqFile}
\alias{readSeqFile}

\alias{readSeqFile.FastqDataFile}
\alias{readFASTQSummary}
\alias{readFASTQSummary.FastqDataFile}

\title{Reads and summarizes a FASTQ file}

\description{
  Reads and summarizes a FASTQ file
}

\usage{
readSeqFile(...)
}

\arguments{
 \item{...}{Arguments passed to \code{\link[qrqc]{readSeqFile}}.}
 \item{seed}{An (optional) \code{\link[base]{integer}} specifying the random seed to be
    set before sampling indices.  The random seed is set to its original
    state when exiting.  If \code{\link[base]{NULL}}, it is not set.}
 \item{cache}{If \code{\link[base:logical]{TRUE}}, memoization is used, otherwise not.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link[qrqc]{FASTQSummary-class}} object.
}

\seealso{
  Internally \code{\link[qrqc]{readSeqFile}} is utilized.
}

\author{Henrik Bengtsson}


\keyword{methods}
\keyword{internal}
