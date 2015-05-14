%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findExternal.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findExternal}
\alias{findExternal}

\alias{getExternalHome}
\alias{findFastqDump}
\alias{findFastQC}
\alias{findJava}
\alias{findPerl}
\alias{findPicard}
\alias{findPython}
\alias{findBowtie2}
\alias{findBWA}
\alias{findCNVkit}
\alias{findGATK}
\alias{findHTSeq}
\alias{findSamtools}
\alias{findTopHat}
\alias{findTopHat1}
\alias{findTopHat2}

\title{Locates an external executable}

\description{
 Locates an external executable.
}

\usage{
  # The generic internal function used
  findExternal(mustExist=TRUE, command, path=NULL, version=NULL, versionPattern=NULL,
  expectedStatus=c(0L, 1L), force=FALSE, verbose=FALSE, ...)

  # Programming environments
  findJava(...)
  findPerl(...)
  findPython(...)

  # fastq-dump
  findFastqDump(..., commandName="fastq-dump",
  versionPattern=".*fastq-dump[ ]*:[ ]*([0-9.]+)", verbose=FALSE)

  # FastQC
  findFastQC(mustExist=TRUE, ..., verbose=FALSE)

  # BWA
  findBWA(...)

  # Bowtie and TopHat
  findBowtie2(..., command=c("bowtie2", "bowtie2-align", "bowtie2-build",
  "bowtie2-inspect"))
  findTopHat(..., command="tophat", path=NULL)
  findTopHat1(..., command="tophat", version=c(1, 2))
  findTopHat2(..., command="tophat2", version=c(2, 3))

  # GATK
  findGATK(mustExist=TRUE, ..., verbose=FALSE)

  # Picard
  findPicard(mustExist=TRUE, ..., verbose=FALSE)

  # Samtools
  findSamtools(...)

  # CNVkit
  findCNVkit(..., command="cnvkit.py")

  # HTSeq
  findHTSeq(..., command=c("htseq-count", "htseq-qa"))
}

\arguments{
  \item{mustExist}{If \code{\link[base:logical]{TRUE}}, an exception is thrown if the executable
     could not be located.}
  \item{command}{A \code{\link[base]{character}} string specifying the name of the
     executable to locate.}
  \item{version}{If non-\code{\link[base]{NULL}}, specifies which version of the
     executable to retrieve.}
  \item{versionPattern}{(A named \code{\link[base]{character}} string specifying the
     \code{\link[base]{gsub}}() regular expression to extraction the version
     where there name is the command-line option specifying how
     to call the external for retrieving the version output.}
  \item{expectedStatus}{An \code{\link[base]{integer}} \code{\link[base]{vector}} of expected status codes
     returned when querying the executable for the version.}
  \item{force}{If \code{\link[base:logical]{TRUE}}, cached results are ignored, otherwise not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{...}{Additional arguments passed to \code{\link{findExternal}}(), or ignored.}
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
