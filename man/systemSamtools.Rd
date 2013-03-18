%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemSamtools.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{systemSamtools}
\alias{systemSamtools.default}
\alias{systemSamtools}

 \title{Calls the samtools executable}

 \description{
  Calls the samtools executable.
 }

 \usage{
\method{systemSamtools}{default}(command, ..., .fake=FALSE, verbose=FALSE)
}

 \arguments{
   \item{command}{A \code{\link[base]{character}} string specifying the samtools command.}
   \item{...}{Additional arguments specifying samtools command line switches.}
   \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \author{Henrik Bengtsson}

