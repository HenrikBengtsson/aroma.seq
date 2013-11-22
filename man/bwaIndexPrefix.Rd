%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  bwaIndexPrefix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{bwaIndexPrefix}
\alias{bwaIndexPrefix.default}
\alias{bwaIndexPrefix}

\title{Generates a prefix for the index files}

\description{
 Generates a prefix for the index files based on a FASTA pathname.
}

\usage{
\method{bwaIndexPrefix}{default}(pathnameFA, subdir="bwa", tags="*", ...)
}

\arguments{
  \item{pathnameFA}{The FASTA file.}
  \item{subdir}{The subdirectory relative to the FASTA file where to put
    the BWA index files.}
  \item{tags}{Tags added to the directory of the index set.}
  \item{...}{Not used.}
}

\examples{
  pathnameFA <- "annotationData/organisms/LambdaPhage/lambda_virus.fa"
  prefix <- bwaIndexPrefix(pathnameFA)
  print(prefix)
}

\author{Henrik Bengtsson}


\keyword{internal}
