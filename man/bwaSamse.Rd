%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  bwaSamse.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{bwaSamse}
\alias{bwaSamse.default}
\alias{bwaSamse}

 \title{Calls the BWA 'samse' command}

 \description{
  Calls the BWA 'samse' command.
 }

 \usage{\method{bwaSamse}{default}(pathnameSAI, pathnameFQ, indexPrefix, pathnameD, ..., verbose=FALSE)}

 \arguments{
   \item{pathnameSAI}{The SAI file to be aligned.}
   \item{pathnameFQ}{The FASTQ file to be aligned.}
   \item{indexPrefix}{The pathname prefix to the BWA index files.}
   \item{pathnameD}{The destination pathname.}
   \item{...}{Additional arguments specifying BWA 'samse' switches
     passed to \code{\link{systemBWA}}().}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \examples{\dontrun{
   pathnameFA <- "annotationData/organisms/LambdaPhage/lambda_virus.fa"
   bwaIndex(pathnameFA)

   pathnameSAI <- "bwaData/LambdaVirusExample/Generic/reads_1.sai";
   pathnameFQ <- "fastqData/LambdaVirusExample/Generic/reads_1.fq";
   pathnameD <- "bwaData/LambdaVirusExample/Generic/reads_1.sam";
   bwaSamse(pathnameSAI=pathnameSAI, pathnameFQ=pathnameFQ, pathnameFA=pathnameFA, pathnameD=pathnameD);
 }}

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

