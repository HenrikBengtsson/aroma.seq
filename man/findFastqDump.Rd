%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findFastqDump.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findFastqDump}
\alias{findFastqDump}


\title{Locates the fastq-dump executable}

\description{
 Locates the fastq-dump executable on the current system.
}

\usage{
findFastqDump(..., commandName="fastq-dump",
  versionPattern=".*fastq-dump[ ]*:[ ]*([0-9.]+)", verbose=FALSE)
}

\arguments{
  \item{...}{Not used.}
  \item{commandName}{'fastq-dump'; command name to find}
  \item{versionPattern}{regexp to use if version not found properly; default should work}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{internal}
