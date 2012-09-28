%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemGATK.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{systemGATK}
\alias{systemGATK.default}
\alias{systemGATK}

 \title{Calls the GATK executable}

 \description{
  Calls the GATK executable.
 }

 \usage{\method{systemGATK}{default}(..., verbose=FALSE)}

 \arguments{
   \item{...}{Arguments specifying GATK command line switches.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \examples{\dontrun{
   pathnameBAM <- "bwaData/LambdaVirusExample,bwa,is/Generic/reads_1.bam"
   res <- systemGATK(T="CountReads", ..., stderr=FALSE)
 }}


 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

