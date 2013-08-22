%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doBowtie2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{doBowtie2}
\alias{doBowtie2.default}
\alias{doBowtie2}
\alias{doBowtie2.FastqDataSet}

\title{Read alignment using the Bowtie v2 aligner}

\description{
 Read alignment using the Bowtie v2 aligner based on [1].
}

\usage{
  doBowtie2(...)
  \method{doBowtie2}{FastqDataSet}(dataSet, reference, ..., verbose=FALSE)
}

\arguments{
 \item{dataSet, df}{A \code{\link{FastqDataSet}}.}
 \item{reference}{A \code{\link{FastaReferenceFile}} or a \code{\link{Bowtie2IndexSet}} specifying the genome reference to align the FASTQ reads to.}
 \item{...}{Additional arguments passed to \code{\link{Bowtie2Alignment}}.}
 \item{verbose}{See \code{\link{Verbose}}.}
}

\value{
  Returns a \code{\link{BamDataSet}}.
}

\references{
 [1] Bowtie2, John Hopkins University, 2013.
     \url{http://bowtie-bio.sourceforge.net/bowtie2/}
}

\author{Henrik Bengtsson}

\seealso{
 For more details, see \code{\link{Bowtie2Alignment}}.
}


\keyword{internal}
