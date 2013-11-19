%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doQDNAseq.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{doQDNAseq}
\alias{doQDNAseq.default}
\alias{doQDNAseq}
\alias{doQDNAseq.BamDataFile}
\alias{doQDNAseq.BamDataSet}
\alias{doQDNAseq.FastqDataSet}

\title{Quantitative inference of copy number aberrations with DNA isolated from fresh or formalin-fixed tissues by shallow whole-genome sequencing (QDNAseq)}

\description{
 Quantitative inference of copy number aberrations with DNA isolated from fresh or formalin-fixed tissues by shallow whole-genome sequencing (QDNAseq) based on [1].
 The algorithm is processed in bounded memory, meaning virtually
 any number of samples can be analyzed on also very limited computer
 systems.
}

\usage{
  \method{doQDNAseq}{FastqDataSet}(dataSet, binWidth, reference, ..., verbose=FALSE)
  \method{doQDNAseq}{BamDataSet}(dataSet, ..., force=FALSE, verbose=FALSE)
  \method{doQDNAseq}{BamDataFile}(df, binWidth, residual=TRUE, blacklist=TRUE, mappability=NA, bases=NA,
  filterAllosomes=TRUE, ..., path=".", force=FALSE, verbose=FALSE)
}

\arguments{
 \item{dataSet, df}{A \code{\link{FastqDataSet}} or a \code{\link{BamDataSet}} (or a \code{\link{BamDataFile}}.}
 \item{binWidth}{A positive \code{\link[base]{numeric}} specifying the bin width (in units of kbp).
   Alternatively, a \code{\link[Biobase]{AnnotatedDataFrame}} specifying the bins.}
 \item{reference}{A \code{\link{FastaReferenceFile}} or a \code{\link{BwaIndexSet}} specifying the genome reference to align the FASTQ reads to.}
 \item{...}{Additional arguments passed to \code{\link[QDNAseq]{applyFilters}},
   \code{\link[QDNAseq]{correctBins}} and \code{\link[QDNAseq]{normalizeBins}}.}
 \item{force}{If \code{\link[base:logical]{TRUE}}, cached results are ignored.}
 \item{verbose}{See \code{\link{Verbose}}.}
}

\value{
  Returns a \code{\link[R.filesets]{GenericDataFileSet}} containing
  \code{\link[QDNAseq]{QDNAseqReadCounts}} objects.
}

\references{
 [1] TBA.
}

\author{Henrik Bengtsson}


\keyword{internal}
