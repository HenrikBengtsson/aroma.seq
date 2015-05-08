%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  createIndexPrefix.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{createIndexPrefix}
\alias{createIndexPrefix.default}
\alias{createIndexPrefix}

\title{Generates a prefix for a set of index files}

\description{
 Generates a prefix for a set of index files based on a FASTA pathname.
}

\usage{
\method{createIndexPrefix}{default}(pathnameFA, subdir=NULL, tags="*", asteriskTags=NULL, ...)
}

\arguments{
  \item{pathnameFA}{The FASTA file.}
  \item{subdir}{The subdirectory relative to the FASTA file where to put
    the index files.}
  \item{tags}{Tags added to the directory of the index set.}
  \item{asteriskTags}{Tags to replace \code{"*"} in argument \code{tags}.}
  \item{...}{Not used.}
}

\examples{
  pathnameFA <- "annotationData/organisms/Lambda_phage/lambda_virus.fa"
  prefix <- createIndexPrefix(pathnameFA)
  print(prefix)
  prefix <- createIndexPrefix(pathnameFA, tags="*,foo")
  print(prefix)
}

\author{Henrik Bengtsson}


\keyword{internal}
