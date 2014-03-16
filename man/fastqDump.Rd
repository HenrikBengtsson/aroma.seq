%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  fastqDump.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{fastqDump}
\alias{fastqDump.default}
\alias{fastqDump}

\title{Calls the fastq-dump executable}

\description{
 Calls the fastq-dump executable.
}

\usage{
\method{fastqDump}{default}(..., split3=TRUE, outPath="fastqData", pathnames=character(0L), verbose=FALSE)
}

\arguments{
  \item{...}{Additional arguments passed to \code{\link{systemFastqDump}}().}
  \item{outPath}{Directory where output fastq files will reside at completion}
  \item{pathnames}{Zero or more input pathnames.}
  \item{split3}{If TRUE, add '--split-3' option to produce two files in case of paired-end data}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{internal}
