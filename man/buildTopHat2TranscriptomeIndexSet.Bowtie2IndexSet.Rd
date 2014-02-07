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

\title{Calls the TopHat executable to build a transcriptome index on a gene model in a reference genome}

\description{
 Calls the TopHat executable to build a transcriptome index on a gene model in a reference genome.
}

\arguments{
  \item{tOutDir}{(optional) Output dir for tophat run itself (log files)}
  \item{gtf}{GtfDataFile to be indexed}
  \item{tiOutPrefix}{(optional) Output path and prefix for transcriptome index}
  \item{...}{(Not used)}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}

\references{
 [1] TopHat, University of Maryland, 2013.
     \url{http://http://tophat.cbcb.umd.edu/}
}


\keyword{internal}
\keyword{methods}
