%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  bowtie2.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{bowtie2}
\alias{bowtie2}


\title{Calls the Bowtie2 executable to align input reads}

\description{
 Calls the Bowtie2 executable to align input reads.
}

\usage{
bowtie2(reads1, reads2=NULL, indexPrefix, pathnameSAM, ..., gzAllowed=NA, verbose=FALSE)
}

\arguments{
  \item{reads1}{(required) A \code{\link[base]{vector}} of FASTQ pathnames of reads.}
  \item{reads2}{(optional; paired-end only) A \code{\link[base]{vector}} of FASTQ pathnames of mate reads.}
  \item{indexPrefix}{Bowtie2 reference index prefix.}
  \item{pathnameSAM}{Output SAM file.}
  \item{...}{...}
  \item{gzAllowed}{A \code{\link[base]{logical}} specifying whether gzipped FASTQ files are supported or not.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns ...
}

\author{Henrik Bengtsson}


\keyword{internal}
