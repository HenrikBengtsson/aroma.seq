%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findBWA.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{findBWA}
\alias{findBWA.default}
\alias{findBWA}

 \title{Locates the BWA executable}

 \description{
  Locates the BWA executable on the current system.
 }

 \usage{\method{findBWA}{default}(mustExists=TRUE, ..., verbose=FALSE)}

 \arguments{
   \item{mustExists}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
      could not be located.}
   \item{...}{Not used.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \details{
  The BWA executable is searched for as follows:
  \enumerate{
   \item \code{Sys.which("bwa")}
  }
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

