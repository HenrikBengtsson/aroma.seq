%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fastQC.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fastQC}
\alias{fastQC.default}
\alias{fastQC}

\title{Calls the FastQC executable}

\description{
 Calls the FastQC executable.
}

\usage{
\method{fastQC}{default}(..., pathnames=character(0L), outPath=NULL, verbose=FALSE)
}

\arguments{
  \item{pathnames}{Zero or more input pathnames.}
  \item{...}{Additional arguments passed to \code{\link{systemFastQC}}().}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Henrik Bengtsson}

\references{
 [1] Simon Andrews,
     FastQC - A quality control tool for high throughput sequence data,
     March 2014.
     \url{http://www.bioinformatics.babraham.ac.uk/projects/fastqc/}
 [2] Simon Andrews,
     FastQC: A quality control application for FastQ data [ANNOUNCEMENT],
     2010-04-26.
     \url{http://seqanswers.com/forums/showthread.php?t=4846}
}


\keyword{internal}
