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

\usage{\method{systemCuffdiff}{default}(commandName="cuffdiff", gtfFile, samFiles, ..., system2ArgsList=list(stdout = FALSE), .fake=FALSE, verbose=FALSE)}

\arguments{
  \item{...}{Additional arguments specifying Cuffdiff command line switches.}
  \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{TBA}

