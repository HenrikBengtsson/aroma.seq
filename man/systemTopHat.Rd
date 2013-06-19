%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemTopHat.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{systemTopHat}
\alias{systemTopHat.default}
\alias{systemTopHat}

\title{Wrapper for tophat}

\description{
 Wrapper for tophat.
}

\usage{
\method{systemTopHat}{default}(commandName="tophat", ...,
  system2ArgsList=list(stdout = TRUE, stderr = FALSE), .fake=FALSE, verbose=FALSE)
}

\arguments{
  \item{commandName}{A \code{\link[base]{character}} string specifying the tophat command.}
  \item{...}{Additional arguments specifying TopHat command line switches.}
  \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
  \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{internal}
