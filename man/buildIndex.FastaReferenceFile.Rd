%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastaReferenceFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{buildIndex.FastaReferenceFile}
\alias{buildIndex.FastaReferenceFile}
\alias{FastaReferenceFile.buildIndex}
\alias{buildIndex.FastaReferenceFile}
\alias{buildIndex,FastaReferenceFile-method}


 \title{Builds an FAI index file}

 \description{
   Builds an FAI index file.
 }

 \usage{\method{buildIndex}{FastaReferenceFile}(this, ..., skip=TRUE, verbose=FALSE)}

 \arguments{
  \item{...}{Additional arguments passed to \code{\link[Rsamtools]{indexFa}}.}
  \item{skip}{If \code{\link[base:logical]{TRUE}}, the index files are not rebuilt if already available.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns a \code{\link[R.filesets]{GenericDataFileSet}} consisting of the BWA index files.
 }

 \author{TBA}
\keyword{internal}
\keyword{methods}
