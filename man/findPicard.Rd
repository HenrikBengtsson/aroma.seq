%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findPicard.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findPicard}
\alias{findPicard}

\alias{PICARD_HOME}

\title{Locates the Picard executable}

\description{
 Locates the Picard executable on the current system.
}

\usage{
findPicard(mustExists=TRUE, ..., verbose=FALSE)
}

\arguments{
  \item{mustExists}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\details{
 The Picard tools directory is searched for as follows:
 \enumerate{
  \item \code{Sys.getenv("PICARD_HOME")}
 }
}

\author{Henrik Bengtsson}

