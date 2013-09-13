%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  bowtie2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{bowtie2}
\alias{bowtie2.default}
\alias{bowtie2}
\alias{bowtie2_hb}

\title{Calls the Bowtie2 executable to align input reads}

\description{
 Calls the Bowtie2 executable to align input reads.
}

\usage{
\method{bowtie2}{default}(bowtieRefIndexPrefix=NULL, reads1=NULL, reads2=NULL, readsU=NULL, samFile=NULL,
  optionsVec, ..., commandName="bowtie2", verbose=FALSE)
}

\arguments{
  \item{bowtieRefIndexPrefix}{bowtie2 reference index (partial pathname, minus the .1.bt2 suffix)}
  \item{reads1}{Vector of fastq files to align, paired with reads2}
  \item{reads2}{Vector of fastq files to align, paired with reads1}
  \item{readsU}{Vector of fastq files to align (at least one of reads1 or readsU must be non-null}
  \item{samFile}{Output file name}
  \item{optionsVec}{Vector of named options (do not include names x, 1, 2, U, or S)}
  \item{...}{...}
  \item{commandName}{Name of executable}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}


\author{Taku Tokuyasu, Henrik Bengtsson}


\keyword{internal}
