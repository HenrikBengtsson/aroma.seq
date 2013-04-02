%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findTopHat.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findTopHat}
\alias{findTopHat.default}
\alias{findTopHat}
\alias{findTopHat1}
\alias{findTopHat1.default}
\alias{findTopHat2}
\alias{findTopHat2.default}

\title{Locates the TopHat executable of a certain version}

\description{
 Locates the TopHat executable of a certain version on the current system.
}

\usage{
\method{findTopHat}{default}(mustExists=TRUE, command="tophat", version, ..., verbose=FALSE)
}

\arguments{
  \item{mustExists}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{command}{A \code{\link[base]{character}} string specifying the name of the executable.}
  \item{version}{If non-\code{\link[base]{NULL}}, specifies which version of TopHat.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\details{
 The TopHat executable is searched for as follows:
 \enumerate{
  \item \code{Sys.which("tophat")}
 }
}

\author{Henrik Bengtsson}

