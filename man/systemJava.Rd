%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemJava.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{systemJava}
\alias{systemJava.default}
\alias{systemJava}

\title{Calls the Java executable}

\description{
 Calls the Java executable.
}

\usage{
\method{systemJava}{default}(..., .fmtArg=c(`( |=)$` = "\%s\%s", `.*` = "\%s \%s"), .fake=FALSE, verbose=FALSE)
}

\arguments{
  \item{...}{Additional arguments specifying Java command line switches.}
  \item{.fmtArg}{A named \code{\link[base]{character}} \code{\link[base]{vector}} specifying how command line
   options are formatted, e.g. \code{"-f bar"} vs  \code{"--foo=bar"}.}
  \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Henrik Bengtsson}


\keyword{internal}
