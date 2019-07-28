%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  gatk.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{gatk}
\alias{gatk}

\alias{gatk.default}

\title{Calls the GATK executable}

\description{
 Calls the GATK executable.
}

\usage{
 \method{gatk}{default}(analysisType, ..., pathnameI=NULL, pathnameR=NULL, pathnameL=NULL, pathnameLog=NULL,
  outPath=NULL, overwrite=FALSE, verbose=FALSE)
}

\arguments{
  \item{analysisType}{(required) A \code{\link[base]{character}} string specifying type of
    analysis.
    In GATK, this corresponds to engine argument
    \code{--analysis_type} (or \code{-T}).}
  \item{...}{(optional) named arguments}.
  \item{pathnameI}{(optional) A \code{\link[base]{character}} \code{\link[base]{vector}} specifying one or more
    (SAM or BAM) files.
    In GATK, this corresponds to one or more engine arguments
    \code{--input_file} (or \code{-I}).}
  \item{pathnameR}{(optional) A \code{\link[base]{character}} string specifying a reference file.
    In GATK, this corresponds to engine argument
    \code{--reference_sequence} (or \code{-R}).}
  \item{pathnameL}{(optional) A \code{\link[base]{character}} string specifying a ROD file.
    In GATK, this corresponds to engine argument
    \code{--interval} (or \code{-L}).}
  \item{outPath}{Directory where result files are written.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Henrik Bengtsson}

\references{
 [1] GATK: The Genome Analysis Toolkit,
     Broad Institute, 2014.
     \url{http://www.broadinstitute.org/gatk/}
}


\keyword{methods}
\keyword{internal}
