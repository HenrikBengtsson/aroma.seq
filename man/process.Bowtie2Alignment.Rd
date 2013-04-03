%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Bowtie2Alignment.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{process.Bowtie2Alignment}
\alias{process.Bowtie2Alignment}
\alias{Bowtie2Alignment.process}
\alias{process,Bowtie2Alignment-method}

 \title{Runs the aligner}

 \description{
   Runs the aligner on all input files.
   The generated BAM files are all sorted and indexed.
 }

 \usage{
\method{process}{Bowtie2Alignment}(this, ..., skip=TRUE, force=FALSE, verbose=FALSE)
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
