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
\method{systemCufflinks}{default}(commandName="cufflinks", ..., system2ArgsList=list(stdout = FALSE), .fake=FALSE, verbose=FALSE)
}

\arguments{
  \item{...}{Additional arguments specifying Cufflinks command line switches.}
  \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Henrik Bengtsson}

