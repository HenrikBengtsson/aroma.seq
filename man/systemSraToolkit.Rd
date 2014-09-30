%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemSraToolkit.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{systemSraToolkit}
\alias{systemSraToolkit.default}
\alias{systemSraToolkit}

\title{Calls an SRA Toolkit executable}

\description{
 Calls an SRA Toolkit executable.
}

\usage{
\method{systemSraToolkit}{default}(command, ..., verbose=FALSE)
}

\arguments{
  \item{command}{The SRA Toolkit executable to call.}
  \item{...}{Command line switches.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\examples{\dontrun{
  pathnameSRA <- "sraData/DataSetA/HomoSapiens/reads.sra"
  res <- systemSraToolkit("fastq-dump", pathnameSRA, ..., stderr=FALSE)
}}


\author{Henrik Bengtsson}


\keyword{internal}
