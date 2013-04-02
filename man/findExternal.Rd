%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findExternal.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findExternal}
\alias{findExternal}

\alias{findJava}
\alias{findPerl}
\alias{findBowtie2}
\alias{findBWA}
\alias{findSamtools}
\alias{findTopHat}
\alias{findTopHat1}
\alias{findTopHat2}

\title{Locates an external executable}

\description{
 Locates an external executable.
}

\usage{
  findJava(...)
  findPerl(...)
  findBowtie2(..., command=c("bowtie2", "bowtie2-align", "bowtie2-build", "bowtie2-inspect"))
  findBWA(...)
  findSamtools(...)
  findTopHat(..., command=tophat)
  findTopHat1(..., command=tophat, version=c(1, 2))
  findTopHat2(..., command=tophat2, version=c(2, 3))
}

\arguments{
  \item{mustExist}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{command}{A \code{\link[base]{character}} string specifying the name of the
     executable to locate.}
  \item{version}{If non-\code{\link[base]{NULL}}, specifies which version of the
     executable to trieve.}
  \item{versionPattern}{(A named \code{\link[base]{character}} string specifying the
     \code{\link[base]{gsub}}() regular expression to extraction the version
     where there name is the command-line option specifying how
     to call the external for retrieving the version output.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, cached results are ignored, otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{...}{Not used.}
}

\value{
  Returns the pathname (or the path) of the external executable.
  If not found, \code{\link[base]{NULL}} is returned, unless if \code{mustExist=TRUE}
  in case an error is thrown.
  If \code{versionPattern} is specified, then the inferred version
  is returned as attribute 'version'.
}

\details{
 The executable is searched using (in order):
 \enumerate{
  \item \code{Sys.which(command)}
 }
}

\author{Henrik Bengtsson}


\keyword{internal}
