%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemCuffdiff.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{systemCuffdiff}
\alias{systemCuffdiff.default}
\alias{systemCuffdiff}

\title{Calls the Cuffdiff executable to perform isoform abundance estimation}

\description{
 Calls the Cuffdiff executable to perform isoform abundance estimation.
}

\usage{
\method{systemCuffdiff}{default}(cuffdiffStr="", system2ArgsList=list(stdout = FALSE), ...,
  commandName="cuffdiff", .fake=TRUE, verbose=TRUE)
}

\arguments{
  \item{cuffdiffStr}{A \code{\link[base]{character}} string of arguments to the cuffdiff command}
  \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
  \item{...}{...}
  \item{commandName}{Name of executable}
  \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{internal}
