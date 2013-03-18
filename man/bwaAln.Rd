%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  bwaAln.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{bwaAln}
\alias{bwaAln.default}
\alias{bwaAln}

 \title{Calls the BWA 'aln' command}

 \description{
  Calls the BWA 'aln' command.
 }

 \usage{
\method{bwaAln}{default}(pathnameFQ, indexPrefix, pathnameD, ..., verbose=FALSE)
}

 \arguments{
   \item{pathnameFQ}{The FASTQ file to be aligned.}
   \item{indexPrefix}{The pathname prefix to the BWA index files.}
   \item{...}{Additional arguments specifying BWA 'aln' switches
     passed to \code{\link{systemBWA}}().}
   \item{pathnameD}{The destination pathname.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \examples{\dontrun{
   pathnameFA <- "annotationData/organisms/LambdaPhage/lambda_virus.fa"
   bwaIndex(pathnameFA, method="is")
   indexPrefix <- bwaIndexPrefix(pathnameFA, method="is")
   bwaAln("fastqData/LambdaVirusExample/Generic/reads_1.fq", indexPrefix=indexPrefix, pathnameD="fastqData/LambdaVirusExample/Generic/reads_1.sai")
 }}

 \author{Henrik Bengtsson}

