%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findBowtie2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findBowtie2}
\alias{findBowtie2.default}
\alias{findBowtie2}

\title{Locates one of the bowtie2 executable}

\description{
 Locates one of the bowtie2 executable on the current system.
}

\usage{
\method{findBowtie2}{default}(mustExists=TRUE, ..., command=c("bowtie2", "bowtie2-align", "bowtie2-build", "bowtie2-inspect"), verbose=FALSE)
}

\arguments{
  \item{mustExists}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{...}{Not used.}
  \item{command}{A \code{\link[base]{character}} string specifying the executable to locate.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\details{
 The Bowtie2 executable is searched for as follows:
 \enumerate{
  \item \code{Sys.which(command)}
 }
}

\author{Henrik Bengtsson}

