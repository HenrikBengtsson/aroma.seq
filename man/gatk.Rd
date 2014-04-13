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
 \method{gatk}{default}(..., pathnameI=NULL, pathnameR=NULL, analysisType, outPath="gatkData/",
  verbose=FALSE)
}

\arguments{
  \item{...}{(optional) named arguments}.
  \item{analysisType}{(required) A \code{\link[base]{character}} string specifying type of analysis.
    In GATK, this corrensponds to argument \code{--analysis_type} (or \code{-T}).}
  \item{pathnameI}{(optional) A \code{\link[base]{character}} \code{\link[base]{vector}} specifying one or more
    (SAM or BAM) files.
    In GATK, this corrensponds to one or more arguments
    \code{--input_file} (or \code{-I}).}
  \item{pathnameR}{(optional) A \code{\link[base]{character}} string specifying an reference file.
    In GATK, this corrensponds to argument \code{--reference_sequence} (or \code{-R}).}
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
