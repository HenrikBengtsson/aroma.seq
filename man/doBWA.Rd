%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doBWA.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{doBWA}
\alias{doBWA.default}
\alias{doBWA}
\alias{doBWA.FastqDataSet}

\title{Read alignment using the Burrows-Wheeler Transform aligner (BWA)}

\description{
 Read alignment using the Burrows-Wheeler Transform aligner (BWA) based on [1].
}

\usage{
  doBWA(...)
  \method{doBWA}{FastqDataSet}(dataSet, reference, ..., verbose=FALSE)
}

\arguments{
 \item{dataSet, df}{A \code{\link{FastqDataSet}}.}
 \item{reference}{A \code{\link{FastaReferenceFile}} or a \code{\link{BwaIndexSet}} specifying the genome reference to align the FASTQ reads to.}
 \item{...}{Additional arguments passed to \code{\link{BwaAlignment}}.}
 \item{verbose}{See \code{\link{Verbose}}.}
}

\value{
  Returns a \code{\link{BamDataSet}}.
}

\references{
  [1] Li H. and Durbin R., \emph{Fast and accurate short read alignment
      with Burrows-Wheeler Transform}. Bioinformatics, 2009.\cr
  [2] Li H. and Durbin R., \emph{Fast and accurate long-read alignment
      with Burrows-Wheeler Transform}. Bioinformatics, 2010.\cr
}

\author{Henrik Bengtsson}

\seealso{
 For more details, see \code{\link{BwaAlignment}}.
}


\keyword{internal}
