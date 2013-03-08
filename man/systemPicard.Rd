%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemPicard.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{systemPicard}
\alias{systemPicard.default}
\alias{systemPicard}

 \title{Calls the Picard executable}

 \description{
  Calls the Picard executable.
 }

 \usage{\method{systemPicard}{default}(command, ..., verbose=FALSE)}

 \arguments{
   \item{command}{A \code{\link[base]{character}} string specifying the Picard command
     (the name of the jar file without the *.jar extension).}
   \item{...}{Additional arguments specifying Picard command line switches.}
   \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \examples{\dontrun{
   pathnameBAM <- "bwaData/LambdaVirusExample,bwa,is/Generic/reads_1.bam"
   res <- systemPicard("ValidateSamFile", INPUT=pathnameBAM, stderr=FALSE)
   ## ERROR: Read groups is empty
   print(res)
   ## [1] 1

   res <- systemPicard("ValidateSamFile", INPUT=pathnameBAM,
                       IGNORE="MISSING_READ_GROUP", stderr=FALSE)
   ## "No errors found"
   print(res)
   ## [1] 0
 }}


 \author{TBA}

