%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  bwaSampe.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{bwaSampe}
\alias{bwaSampe.default}
\alias{bwaSampe}

\title{Generates BWA-backtrack paired-end (PE) alignments via 'bwa sampe'}

\description{
 Generates BWA-backtrack paired-end (PE) alignments via 'bwa sampe'.
}

\usage{
\method{bwaSampe}{default}(pathnameSAI, pathnameFQ, indexPrefix, pathnameD, ..., verbose=FALSE)
}

\arguments{
  \item{pathnameSAI}{A \code{\link[base]{character}} \code{\link[base]{vector}} of two SAI files.}
  \item{pathnameFQ}{A \code{\link[base]{character}} \code{\link[base]{vector}} of two FASTQ or BAM files.}
  \item{indexPrefix}{The pathname prefix to the BWA index files.}
  \item{pathnameD}{The destination pathname.}
  \item{...}{Additional arguments specifying BWA 'sampe' switches
    passed to \code{\link{systemBWA}}().}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\examples{\dontrun{
  pathnameFA <- "annotationData/organisms/LambdaPhage/lambda_virus.fa"
  bwaIndex(pathnameFA)

  pathnameSAI <- "bwaData/LambdaVirusExample/Generic/reads_1.sai";
  pathnameFQ <- "fastqData/LambdaVirusExample/Generic/reads_1.fq";
  pathnameD <- "bwaData/LambdaVirusExample/Generic/reads_1.sam";
  bwaSampe(pathnameSAI=pathnameSAI, pathnameFQ=pathnameFQ,
           pathnameFA=pathnameFA, pathnameD=pathnameD);
}}

\author{Henrik Bengtsson}


\keyword{internal}
