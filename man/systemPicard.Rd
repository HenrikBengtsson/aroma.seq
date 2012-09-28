%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemPicard.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{systemPicard}
\alias{systemPicard.default}
\alias{systemPicard}

 \title{Calls the Picard executable}

 \description{
  Calls the Picard executable.
 }

 \usage{\method{systemPicard}{default}(command, ..., .fake=FALSE, verbose=FALSE)}

 \arguments{
   \item{command}{A \code{\link[base]{character}} string specifying the Picard command
     (the name of the jar file without the *.jar extension).}
   \item{...}{Additional arguments specifying Picard command line switches.}
   \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

