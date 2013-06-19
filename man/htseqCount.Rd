%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  htseqCount.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{htseqCount}
\alias{htseqCount.default}
\alias{htseqCount}

\title{Calls the htseq-count executable to count input reads on features}

\description{
 Calls the htseq-count executable to count input reads on features.
}

\usage{
\method{htseqCount}{default}(samFile=NULL, gfFile=NULL, outFile=NULL, optionsVec=c(s = "no", a = "10"),
  ..., command="htseq-count", verbose=FALSE)
}

\arguments{
  \item{command}{Name of executable}
  \item{samFile}{sam file containing aligned reads, sorted by *name* if paired end}
  \item{gfFile}{gene feature file, in gff format}
  \item{outFile}{name of file to receive htseq-count output}
  \item{optionsVec}{Vector of named options to pass to htseq-count}
  \item{...}{(Not used)}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}

