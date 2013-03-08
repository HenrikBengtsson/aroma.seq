%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastaReferenceFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{buildBwaIndexSet.FastaReferenceFile}
\alias{buildBwaIndexSet.FastaReferenceFile}
\alias{FastaReferenceFile.buildBwaIndexSet}
\alias{buildBwaIndexSet.FastaReferenceFile}
\alias{buildBwaIndexSet,FastaReferenceFile-method}


 \title{Builds a BWA index files set}

 \description{
   Builds a BWA index files set.
 }

 \usage{\method{buildBwaIndexSet}{FastaReferenceFile}(this, method, ..., skip=TRUE, verbose=FALSE)}

 \arguments{
  \item{method}{A \code{\link[base]{character}} string specifying the algorithm to use.}
  \item{...}{Additional arguments passed to \code{\link{bwaIndex}}().}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the index files are not rebuilt if already available.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns a \code{\link[R.filesets]{GenericDataFileSet}} consisting of the BWA index files.
 }

 \author{TBA}
\keyword{internal}
\keyword{methods}
