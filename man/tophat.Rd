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
 \method{tophat}{default}(bowtieRefIndexPrefix, reads1, reads2=NULL, outPath="tophat/", optionsVec=NULL,
  ..., command="tophat", verbose=FALSE)
 \method{tophat1}{default}(..., command="tophat")
 \method{tophat2}{default}(..., command="tophat2")
}

\arguments{
  \item{bowtieRefIndexPrefix}{bowtie2 reference index (partial pathname, i.e. minus the .x.bt2 suffix)}
  \item{reads1}{(required) Vector of fastq filenames to align; currently only a single filename is supported}
  \item{reads2}{(optional) Vector of fastq filenames to align, paired with reads1; currently only a single filename is supported}
  \item{outPath}{Directory where result files are written}
  \item{optionsVec}{Vector of named options to pass to tophat}
  \item{...}{(Not used)}
  \item{command}{Name of executable}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\section{Support for compressed input files}{
  TopHat (>= 1.3.0) handles FASTQ files that have been compressed by gzip (or bzip2) [1].
  If not supported, this method will give an informative error message about it.
}

\author{Henrik Bengtsson, Taku Tokuyasu}

\references{
 [1] TopHat, University of Maryland, 2013.
     \url{http://http://tophat.cbcb.umd.edu/}
}


\keyword{methods}
\keyword{internal}
