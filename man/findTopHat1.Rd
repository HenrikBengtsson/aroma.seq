%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findTopHat1.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findTopHat1}
\alias{findTopHat1.default}
\alias{findTopHat1}

\title{Locates the TopHat executable and tests version}

\description{
 Locates the TopHat executable and tests version on the current system.
}

\usage{\method{findTopHat1}{default}(mustExists=TRUE, ..., verbose=FALSE)}

\arguments{
  \item{mustExists}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\details{
 The TopHat executable is searched for as follows:
 \enumerate{
  \item \code{Sys.which("tophat")}
 }
}

\author{TBA}

