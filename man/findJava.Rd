%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findJava.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{findJava}
\alias{findJava.default}
\alias{findJava}

 \title{Locates the Java executable}

 \description{
  Locates the Java executable on the current system.
 }

 \usage{\method{findJava}{default}(mustExists=TRUE, ..., verbose=FALSE)}

 \arguments{
   \item{mustExists}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
      could not be located.}
   \item{...}{Not used.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \details{
  The Java executable is searched for as follows:
  \enumerate{
   \item \code{Sys.which("java")}
  }
 }

 \author{TBA}

