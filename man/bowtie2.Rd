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

\title{Calls the Bowtie2 executable to align input reads}

\description{
 Calls the Bowtie2 executable to align input reads.
}

\usage{
\method{bowtie2}{default}(command="bowtie2", bowtieRefIndexPrefix=NULL, reads1=NULL, reads2=NULL, readsU=NULL, samFile=NULL, optionsVec, ..., verbose=FALSE)
}

\arguments{
  \item{command}{Name of executable}
  \item{bowtieRefIndexPrefix}{bowtie2 reference index (partial pathname, i.e. minus the .x.bt2 suffix)}
  \item{reads1}{Vector of fastq files to align, paired with reads2}
  \item{reads2}{Vector of fastq files to align, paired with reads1}
  \item{readsU}{Vector of fastq files to align (at least one of reads1 or readsU must be non-null}
  \item{optionsVec}{Vector of named options (do not include names x, 1, 2, U, or S)}
  \item{...}{...}
}

\examples{\dontrun{
}}

\author{Henrik Bengtsson}

