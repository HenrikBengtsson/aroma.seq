%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findCmdv.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findCmdv}
\alias{findCmdv.default}
\alias{findCmdv}

\title{Locates the executable given by 'command'; tests version if possible}

\description{
 Locates the executable given by 'command'; tests version if possible on the current system.
}

\usage{
\method{findCmdv}{default}(command=NULL, version=NULL, mustExists=TRUE, ..., verbose=FALSE)
}

\arguments{
  \item{command}{Name of executable to find}
  \item{mustExists}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{version}{}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\details{
 The executable is searched for as follows:
 \enumerate{
  \item \code{Sys.which(command)}
 }
 The software version is obtained by trying to parse
 \enumerate{
  \item \code{'cmd --version'}
 }
}

\author{Taku Tokuyasu}


\keyword{internal}
