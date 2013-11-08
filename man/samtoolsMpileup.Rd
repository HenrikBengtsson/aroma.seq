%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  samtoolsMpileup.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{samtoolsMpileup}
\alias{samtoolsMpileup.default}
\alias{samtoolsMpileup}

\title{Calls the samtools 'mpileup' command}

\description{
 Calls the samtools 'mpileup' command.
}

\usage{
\method{samtoolsMpileup}{default}(refFile, bamFile, pathnameD="mpileup.out", ..., verbose=FALSE)
}

\arguments{
  \item{refFile}{input reference file}
  \item{bamFile}{input BAM file (currently only one file is supported)}
  \item{pathnameD}{output file (default 'mpileup.out')}
  \item{...}{Additional arguments specifying samtools 'mpileup' switches
    passed to \code{\link{systemSamtools}}().}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{internal}
