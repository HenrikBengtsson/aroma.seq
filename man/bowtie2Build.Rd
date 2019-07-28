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
\method{bowtie2Build}{default}(pathnameFAs, bowtieRefIndexPrefix, optionsVec=NULL, ...,
  command="bowtie2-build", verbose=FALSE)
}

\arguments{
  \item{pathnameFAs}{A \code{\link[base]{character}} \code{\link[base]{vector}} of FASTA reference files.}
  \item{bowtieRefIndexPrefix}{A \code{\link[base]{character}} string specifying the bowtie2
    reference index to be built (partial pathname, minus the .*.bt2 suffix).}
  \item{optionsVec}{(optional) A named \code{\link[base]{character}} \code{\link[base]{vector}}.}
  \item{...}{...}
  \item{command}{The name of the external executable.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\section{Support for compressed input files}{
  If gzipped FASTA files are used, this function will temporarily decompress
  before passing them to the bowtie2-build external software (which only
  support non-compressed FASTA files).
}

\section{Known issues}{
  The FASTA pathnames must not contain commas.
  If detected, this method generates an informative error.
}

\author{Henrik Bengtsson, Taku Tokuyasu}


\keyword{internal}
