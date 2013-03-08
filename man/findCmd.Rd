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

\title{Locates the executable given by 'command'}

\description{
 Locates the executable given by 'command' on the current system.
}

\usage{\method{findCmd}{default}(command, mustExists=TRUE, ..., verbose=FALSE)}

\arguments{
  \item{command}{Name of executable to find}
  \item{mustExists}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{...}{Not used.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\details{
 The executable is searched for as follows:
 \enumerate{
  \item \code{Sys.which(command)}
 }
 NB: This method does NOT do any version checking!
}

\author{TBA}

