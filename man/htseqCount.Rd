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
\method{htseqCount}{default}(pathnameS, pathnameD=NULL, gff, optionsVec=c(`-s` = "no", `-a` = "10"), ...,
  command="htseq-count", verbose=FALSE)
}

\arguments{
  \item{pathnameS}{source file containing aligned reads; .sam and .bam supported}
  \item{pathnameD}{(Optional) destination file to receive htseq-count output; use pathnameS with .count extension when NULL}
  \item{gff}{gene feature file, in gff format}
  \item{optionsVec}{Vector of named options to htseq-count}
  \item{...}{(Not used)}
  \item{command}{Name of executable}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{internal}
