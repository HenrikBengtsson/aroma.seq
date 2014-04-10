%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doTopHat2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{doTopHat2}
\alias{doTopHat2.default}
\alias{doTopHat2}
\alias{doTopHat2.FastqDataSet}

\title{Read alignment using the TopHat v2 aligner}

\description{
 Read alignment using the TopHat v2 aligner based on [1].
}

\usage{
  doTopHat2(...)
  \method{doTopHat2}{FastqDataSet}(dataSet, reference, transcripts, ..., verbose=FALSE)
}

\arguments{
 \item{dataSet, df}{A \code{\link{FastqDataSet}}.}
 \item{reference}{A \code{\link{FastaReferenceFile}} or a \code{\link{Bowtie2IndexSet}} specifying the genome reference to align the FASTQ reads to.}
 \item{transcripts}{A \code{\link{GtfDataFile}} or \code{\link[base]{NULL}} specifying known transcripts and/or gene model annotations.}
 \item{...}{Additional arguments passed to \code{\link{TopHat2Alignment}}.}
 \item{verbose}{See \code{\link{Verbose}}.}
}

\value{
  Returns a \code{\link{BamDataSet}}.
}

\references{
 [1] TopHat, University of Maryland, 2013.
     \url{http://http://tophat.cbcb.umd.edu/} \cr
 [2] Trapnell et al. \emph{Differential gene and transcript expression
     analysis of RNA-seq experiments with TopHat and Cufflinks}.
     Nat Protoc, 2012.\cr
}

\author{Henrik Bengtsson}

\seealso{
 For more details, see \code{\link{TopHat2Alignment}}.
}


\keyword{internal}
