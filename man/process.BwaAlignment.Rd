%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BwaAlignment.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{process.BwaAlignment}
\alias{process.BwaAlignment}
\alias{BwaAlignment.process}
\alias{process,BwaAlignment-method}

 \title{Runs the BWA aligner}

 \description{
   Runs the BWA aligner on all input files.
   The generated BAM files are all sorted and indexed.
 }

 \usage{
\method{process}{BwaAlignment}(this, ..., skip=TRUE, force=FALSE, verbose=FALSE)
}

 \arguments{
  \item{...}{Not used.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, already processed files are skipped.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns a \code{\link{BamDataSet}}.
 }

 \author{Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
