%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  bowtie2Build.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{bowtie2Build}
\alias{bowtie2Build.default}
\alias{bowtie2Build}

\title{Creates index on reference genome using bowtie2-build}

\description{
 Creates index on reference genome using bowtie2-build.
}

\usage{
\method{bowtie2Build}{default}(refReads=NULL, bowtieRefIndexPrefix=NULL, optionsVec=NULL, ..., commandName="bowtie2-build", verbose=FALSE)
}

\arguments{
  \item{refReads}{Vector of (fasta) files with reference sequence}
  \item{bowtieRefIndexPrefix}{bowtie2 reference index to be built (partial pathname, minus the .*.bt2 suffix)}
  \item{optionsVec}{Vector of named options}
  \item{...}{...}
  \item{commandName}{Name of executable}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}


\author{Henrik Bengtsson}

