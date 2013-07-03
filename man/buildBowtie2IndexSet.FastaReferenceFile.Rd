%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastaReferenceFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{buildBowtie2IndexSet.FastaReferenceFile}
\alias{buildBowtie2IndexSet.FastaReferenceFile}
\alias{FastaReferenceFile.buildBowtie2IndexSet}
\alias{buildBowtie2IndexSet,FastaReferenceFile-method}

 \title{Builds a Bowtie2 index files set}

 \description{
   Builds a Bowtie2 index files set.
 }

 \usage{
\method{buildBowtie2IndexSet}{FastaReferenceFile}(this, ..., skip=TRUE, verbose=FALSE)
}

 \arguments{
  \item{...}{Additional arguments passed to \code{\link{bowtie2Build}}().}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the index files are not rebuilt if already available.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns a \code{\link[R.filesets]{GenericDataFileSet}} consisting of the bowtie2 index files.
 }

% \section{Benchmarking}{
%   Examples of processing times:
%   \itemize{
%    \item human_g1k_v37.fasta: ~120 minutes on System A.
%   }
%   where:
%   \itemize{
%    \item 'System A': one core on a random node on a Linux 64-bit cluster.
%    \item 'System B' is Windows 7 Pro 64-bit on a Lenovo Thinkpad X201.
% }

 \seealso{
   Internally, \code{\link{bowtie2Build}}() is used.
 }

 \author{Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
