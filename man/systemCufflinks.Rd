%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemCufflinks.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{systemCufflinks}
\alias{systemCufflinks.default}
\alias{systemCufflinks}

\title{Calls the Cufflinks executable to perform isoform abundance estimation}

\description{
 Calls the Cufflinks executable to perform isoform abundance estimation.
}

\usage{
\method{systemCufflinks}{default}(commandName="cufflinks", ..., system2ArgsList=list(stdout = FALSE),
  .fake=FALSE, verbose=FALSE)
}

\arguments{
  \item{commandName}{A \code{\link[base]{character}} string specifying the Cufflinks command.}
  \item{...}{Arguments specifying the Cufflinks command line.}
  \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
  \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{internal}
