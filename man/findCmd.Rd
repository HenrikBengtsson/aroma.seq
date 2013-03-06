%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findCmd.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findCmd}
\alias{findCmd.default}
\alias{findCmd}

\title{Locates the executable given by 'cmd'}

\description{
 Locates the executable given by 'cmd' on the current system.
}

\usage{\method{findCmd}{default}(command, mustExists=TRUE, ..., verbose=FALSE)}

\arguments{
  \item{mustExists}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\details{
 The executable is searched for as follows:
 \enumerate{
  \item \code{Sys.which("tophat")}
 }
}

\author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

