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

 \usage{\method{bwaAln}{default}(pathnameFQ, pathnameFA, pathnameD, ..., verbose=FALSE)}

 \arguments{
   \item{pathnameFQ}{The FASTQ file to be aligned.}
   \item{pathnameFA}{The FASTA reference file.}
   \item{...}{Additional arguments specifying BWA 'aln' switches
     passed to \code{\link{systemBWA}}().}
   \item{pathnameD}{The destination pathname.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \examples{\dontrun{
   pathnameFA <- "annotationData/organisms/LambdaPhage/lambda_virus.fa"
   bwaIndex(pathnameFA)
   bwaAln("fastqData/LambdaVirusExample/Generic/reads_1.fq", pathnameD="fastqData/LambdaVirusExample/Generic/reads_1.sai", pathnameFA=pathnameFA)
 }}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

