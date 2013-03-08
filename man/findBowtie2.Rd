%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findBowtie2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{findBowtie2}
\alias{findBowtie2.default}
\alias{findBowtie2}

 \title{Locates the Bowtie2 executable}

 \description{
  Locates the Bowtie2 executable on the current system.
 }

 \usage{\method{findBowtie2}{default}(mustExists=TRUE, ..., verbose=FALSE)}

 \arguments{
   \item{mustExists}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
      could not be located.}
   \item{...}{Not used.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \details{
  The Bowtie2 executable is searched for as follows:
  \enumerate{
   \item \code{Sys.which("bowtie2-align")}
  }
 }

 \author{TBA}

