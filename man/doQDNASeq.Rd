%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doQDNASeq.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{doQDNASeq}
\alias{doQDNASeq.default}
\alias{doQDNASeq}
\alias{doQDNASeq.BamDataFile}
\alias{doQDNASeq.BamDataSet}
\alias{doQDNASeq.FastqDataSet}

\title{Quantitative inference of copy number aberrations with DNA isolated from fresh or formalin-fixed tissues by shallow whole-genome sequencing (QDNASeq)}

\description{
 Quantitative inference of copy number aberrations with DNA isolated from fresh or formalin-fixed tissues by shallow whole-genome sequencing (QDNASeq) based on [1].
 The algorithm is processed in bounded memory, meaning virtually
 any number of samples can be analyzed on also very limited computer
 systems.
}

\usage{
  \method{doQDNASeq}{FastqDataSet}(dataSet, binWidth, reference, ..., verbose=FALSE)
  \method{doQDNASeq}{BamDataSet}(dataSet, binWidth, ..., verbose=FALSE)
  \method{doQDNASeq}{BamDataFile}(df, binWidth, log=TRUE, mappability=50, blacklist=0, residual=2, bases=0, ...,
  force=FALSE, verbose=FALSE)
}

\arguments{
 \item{dataSet, df}{A \code{\link{FastqDataSet}} or a \code{\link{BamDataSet}} (or a \code{\link{BamDataFile}}.}
 \item{binWidth}{A \code{\link[base]{numeric}} specifying the bin width (in units kbp).}
 \item{reference}{A \code{\link{FastaReferenceFile}} or a \code{\link{BwaIndexSet}} specifying the genome reference to align the FASTQ reads to.}
 \item{log}{If \code{\link[base:logical]{TRUE}}, the copy numbers are calculated on the log2 scale.}
 \item{mappability, blacklist, residual, bases}{Post-filter arguments.}
 \item{...}{Ignored, or passed to \code{doQDNASeq()}.}
 \item{force}{If \code{\link[base:logical]{TRUE}}, cached results are ignored.}
 \item{verbose}{See \code{\link{Verbose}}.}
}

\value{
  Returns a \code{\link[qdnaseq]{QDNAseqReadCounts}} object.
}

\references{
 [1] TBA.
}

\author{Henrik Bengtsson}

