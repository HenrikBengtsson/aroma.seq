%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findPerl.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findPerl}
\alias{findPerl.default}
\alias{findPerl}

\title{Locates the Perl executable}

\description{
 Locates the Perl executable on the current system.
}

\usage{
\method{findPerl}{default}(mustExists=TRUE, ..., verbose=FALSE)
}

\arguments{
  \item{mustExists}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\details{
 The Perl executable is searched for as follows:
 \enumerate{
  \item \code{Sys.which("perl")}
 }
}

\author{Henrik Bengtsson}

