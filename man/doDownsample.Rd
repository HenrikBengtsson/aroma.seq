%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doDownsample.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{doDownsample}
\alias{doDownsample.default}
\alias{doDownsample}
\alias{doDownsample.BamDataSet}
\alias{doDownsample.FastqDataSet}

\title{Generates a downsampled FASTQ or BAM data set}

\description{
 Generates a downsampled FASTQ or BAM data set.
}

\usage{
  doDownsample(...)
  \method{doDownsample}{BamDataSet}(dataSet, subset=1e+06, ..., verbose=FALSE)
  \method{doDownsample}{FastqDataSet}(dataSet, subset=subset, ..., verbose=FALSE)
}

\arguments{
 \item{dataSet}{A \code{\link{BamDataSet}} or \code{\link{FastqDataSet}}.}
 \item{subset}{An \code{\link[base]{integer}} specifying the total number of reads to sample,
   or a \code{\link[base]{double}} specifying the fraction of total number of reads to sample.}
 \item{...}{Additional arguments passed to specific downsampler, e.g.
    \code{\link{BamDownsampler}} and \code{\link{FastqDownsampler}}.}
 \item{verbose}{See \code{\link{Verbose}}.}
}

\value{
  Returns a \code{\link{GenericDataFileSet}} of the same class
  as the input data set \code{dataSet}.
}

\author{Henrik Bengtsson}


\keyword{internal}
