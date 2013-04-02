%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findGATK.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findGATK}
\alias{findGATK}

\alias{GATK_HOME}

\title{Locates the GATK executable}

\description{
 Locates the GATK executable on the current system.
}

\usage{
findGATK(mustExist=TRUE, ..., verbose=FALSE)
}

\arguments{
  \item{mustExist}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\details{
 The GATK tools directory is searched for as follows:
 \enumerate{
  \item \code{Sys.getenv("GATK_HOME")}
 }
}

\author{Henrik Bengtsson}


\keyword{internal}
