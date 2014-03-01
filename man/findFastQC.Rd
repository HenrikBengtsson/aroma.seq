%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findFastQC.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findFastQC}
\alias{findFastQC}

\alias{FASTQC_HOME}

\title{Locates the FastQC executable}

\description{
 Locates the FastQC executable on the current system.
}

\usage{
findFastQC(mustExist=TRUE, ..., verbose=FALSE)
}

\arguments{
  \item{mustExist}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\details{
 The FastQC tools directory is searched for as follows:
 \enumerate{
  \item \code{Sys.getenv("FASTQC_HOME")}
 }
}

\author{Henrik Bengtsson}


\keyword{internal}
