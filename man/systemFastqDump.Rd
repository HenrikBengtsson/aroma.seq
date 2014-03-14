%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemFastqDump.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{systemFastqDump}
\alias{systemFastqDump.default}
\alias{systemFastqDump}

\title{Wrapper for fastq-dump SRA utility}

\description{
 Wrapper for fastq-dump SRA utility.
}

\usage{
\method{systemFastqDump}{default}(..., commandName="fastq-dump", Stdout=TRUE, Stderr=FALSE, .fake=FALSE,
  verbose=FALSE)
}

\arguments{
  \item{commandName}{A \code{\link[base]{character}} string specifying the FastqDump command.}
  \item{...}{Additional arguments specifying FastqDump command line switches.}
  \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
  \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{internal}
