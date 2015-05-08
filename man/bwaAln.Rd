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

\title{BWA-backtrack alignment via 'bwa aln'}

\description{
 BWA-backtrack alignment via 'bwa aln'.
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
  pathnameFA <- "annotationData/organisms/Lambda_phage/lambda_virus.fa"
  bwaIndex(pathnameFA)
  indexPrefix <- bwaIndexPrefix(pathnameFA)
  bwaAln("fastqData/LambdaVirusExample/Lambda_phage/reads_1.fq",
         indexPrefix=indexPrefix,
         pathnameD="fastqData/LambdaVirusExample/Lambda_phage/reads_1.sai")
}}

\author{Henrik Bengtsson}


\keyword{internal}
