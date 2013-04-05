%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemGenericCmd.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{systemGenericCmd}
\alias{systemGenericCmd.default}
\alias{systemGenericCmd}

\title{Wrapper to call 'command' using system2()}

\description{
 Wrapper to call 'command' using system2().
}

\usage{
\method{systemGenericCmd}{default}(commandName=NULL, ..., system2ArgsList=list(), .fake=FALSE, verbose=FALSE)
}

\arguments{
  \item{commandName}{Name of command to run}
  \item{...}{Arguments that specify the command line string.}
  \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
  \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{internal}
