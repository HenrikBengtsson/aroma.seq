%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  buildTopHat2TranscriptomeIndexSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{buildTopHat2TranscriptomeIndexSet.Bowtie2IndexSet}
\alias{buildTopHat2TranscriptomeIndexSet.Bowtie2IndexSet}
\alias{Bowtie2IndexSet.buildTopHat2TranscriptomeIndexSet}
\alias{buildTopHat2TranscriptomeIndexSet,Bowtie2IndexSet-method}
\alias{buildTopHat2TranscriptomeIndexSet}

\title{Calls TopHat to build a transcriptome index}

\description{
 Calls TopHat to build a transcriptome index.
}

\arguments{
  \item{gtf}{GtfDataFile to be indexed}
  \item{outPath}{(optional) Output directory for index and log file.}
  \item{tiPrefix}{(optional) Prefix for transcriptome index.}
  \item{...}{Arguments passed to tophat().}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the index files are not rebuilt if already available.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}

\references{
 [1] TopHat, University of Maryland, 2013.
     \url{http://http://tophat.cbcb.umd.edu/}
}


\keyword{internal}
\keyword{methods}
