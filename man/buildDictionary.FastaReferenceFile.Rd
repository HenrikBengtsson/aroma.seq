%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastaReferenceFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{buildDictionary.FastaReferenceFile}
\alias{buildDictionary.FastaReferenceFile}
\alias{FastaReferenceFile.buildDictionary}
\alias{buildDictionary,FastaReferenceFile-method}

\title{Builds a DICT sequence dictionary file}

\description{
  Builds a DICT sequence dictionary file.
}

\usage{
\method{buildDictionary}{FastaReferenceFile}(this, ..., skip=TRUE, verbose=FALSE)
}

\arguments{
 \item{...}{Additional arguments passed to \code{\link[Rsamtools]{indexFa}}.}
 \item{skip}{If \code{\link[base:logical]{TRUE}}, the dictionary is not rebuilt if already available.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link[R.filesets]{GenericDataFileSet}} consisting of the BWA index files.
}

\references{
 [1] Geraldine van der Auwera,
     How can I prepare a FASTA file to use as reference?,
     GATK Forum, Sept 2013.
     \url{http://gatkforums.broadinstitute.org/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference}
}

\author{Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
