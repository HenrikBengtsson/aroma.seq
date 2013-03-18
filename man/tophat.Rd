%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  tophat.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{tophat}
\alias{tophat.default}
\alias{tophat}

\title{Calls the TopHat executable to align input reads}

\description{
 Calls the TopHat executable to align input reads.
}

\usage{
\method{tophat}{default}(bowtieRefIndexPrefix=NULL, reads1=NULL, reads2=NULL, optionsVec=NULL, ..., command="tophat", verbose=FALSE)
}

\arguments{
  \item{command}{Name of executable}
  \item{bowtieRefIndexPrefix}{bowtie2 reference index (partial pathname, i.e. minus the .x.bt2 suffix)}
  \item{reads1}{(required) Vector of fastq files to align}
  \item{reads2}{(optional) Vector of fastq files to align, paired with reads1}
  \item{optionsVec}{Vector of named options to pass to tophat}
  \item{...}{(Not used)}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\examples{\dontrun{
}}

\author{Henrik Bengtsson}

