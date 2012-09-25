%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  bwaIndex.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{bwaIndex}
\alias{bwaIndex.default}
\alias{bwaIndex}
 \alias{bwaIndex.FastaReferenceFile}

 \title{Calls the BWA index command}

 \description{
  Calls the BWA index command.
 }

 \usage{\method{bwaIndex}{default}(pathnameFA, indexPrefix="*", ..., verbose=FALSE)}

 \arguments{
   \item{pathnameFA}{The FASTA file to be indexed.}
   \item{indexPrefix}{The prefix for the generated index files.}
   \item{...}{Additional arguments specifying BWA 'index' switches
     passed to \code{\link{systemBWA}}().}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \examples{\dontrun{
   bwaIndex("annotationData/organisms/LambdaPhage/lambda_virus.fa")
 }}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

