%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemBowtie2Build.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{systemBowtie2Build}
\alias{systemBowtie2Build.default}
\alias{systemBowtie2Build}

\title{Wrapper for bowtie2-build}

\description{
 Wrapper for bowtie2-build.
}

\usage{
\method{systemBowtie2Build}{default}(commandName="bowtie2-build", ...,
  system2ArgsList=list(stdout = TRUE, stderr = FALSE), .fake=FALSE, verbose=FALSE)
}

\arguments{
  \item{commandName}{A \code{\link[base]{character}} string specifying the bowtie2 build command.}
  \item{...}{Additional arguments specifying bowtie2-build command line switches.}
  \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
  \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{internal}
