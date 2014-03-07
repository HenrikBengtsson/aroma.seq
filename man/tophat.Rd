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
 \method{tophat}{default}(bowtieRefIndexPrefix, reads1=NULL, reads2=NULL, gtf=NULL, mateInnerDist=NULL,
  mateStdDev=NULL, optionsVec=NULL, ..., outPath="tophat/", command="tophat",
  verbose=FALSE)
 \method{tophat1}{default}(..., command="tophat")
 \method{tophat2}{default}(..., command="tophat2")
}

\arguments{
  \item{bowtieRefIndexPrefix}{A \code{\link[base]{character}} string specifying the Bowtie2 reference index prefix.}
  \item{reads1}{(required) A \code{\link[base]{vector}} of FASTQ pathnames of reads.}
  \item{reads2}{(optional; paired-end only) A \code{\link[base]{vector}} of FASTQ pathnames of mate reads.}
  \item{gtf}{(optional) A GTF pathname.}
  \item{mateInnerDist, mateStdDev}{(optional; paired-end only) The expected mean and standard
        deviation of the inner distance between mate pairs.}
  \item{optionsVec}{Vector of named options to pass to the executable.}
  \item{...}{(Not used)}.
  \item{outPath}{Directory where result files are written.}
  \item{command}{The name of the executable.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\section{Support for compressed input files}{
  TopHat (>= 1.3.0), which was released June 2011, handles FASTQ files that have been compressed
  by gzip (or bzip2) [1]. If not supported, this method will give an informative error message
  about it.
}

\author{Henrik Bengtsson, Taku Tokuyasu}

\references{
 [1] TopHat, University of Maryland, 2013.
     \url{http://http://tophat.cbcb.umd.edu/}
}


\keyword{methods}
\keyword{internal}
