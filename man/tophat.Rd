%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  tophat.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{tophat}
\alias{tophat}

\alias{tophat1}
\alias{tophat2}
\alias{tophat.default}
\alias{tophat1.default}
\alias{tophat2.default}

\title{Calls the TopHat executable to align input reads}

\description{
 Calls the TopHat executable to align input reads.
}

\usage{
 \method{tophat}{default}(bowtieRefIndexPrefix=NULL, reads1=NULL, reads2=NULL, outDir=NULL, optionsVec=NULL,
  ..., .initialTopHatOutDir=NULL, command="tophat", verbose=FALSE)
 \method{tophat1}{default}(..., command="tophat")
 \method{tophat2}{default}(..., command="tophat2")
}

\arguments{
  \item{command}{Name of executable}
  \item{bowtieRefIndexPrefix}{bowtie2 reference index (partial pathname, i.e. minus the .x.bt2 suffix)}
  \item{reads1}{(required) Vector of fastq filenames to align; currently only a single filename is supported}
  \item{reads2}{(optional) Vector of fastq filenames to align, paired with reads1; currently only a single filename is supported}
  \item{.initialTopHatOutDir}{(optional) TopHat output dir used by tophat executable}
  \item{outDir}{(optional) TopHat output dir at method exit}
  \item{optionsVec}{Vector of named options to pass to tophat}
  \item{...}{(Not used)}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{methods}
\keyword{internal}
