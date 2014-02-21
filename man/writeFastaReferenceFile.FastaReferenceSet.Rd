%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastaReferenceSet.writeFastaReferenceFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\name{writeFastaReferenceFile.FastaReferenceSet}
\alias{writeFastaReferenceFile.FastaReferenceSet}
\alias{FastaReferenceSet.writeFastaReferenceFile}
\alias{writeFastaReferenceFile,FastaReferenceSet-method}
\alias{writeFastaReferenceFile}

\title{Writes the content of multiple FASTA files into a new FASTA file}

\description{
  Writes the content of multiple FASTA files into a new FASTA file.  Existing files are left unmodified.
}

\usage{
\method{writeFastaReferenceFile}{FastaReferenceSet}(this, filename=sprintf("\%s,\%s.fa", getOrganism(this),
  getChecksum(this)), path=file.path("annotationData", "organisms", getOrganism(this)),
  ..., overwrite=FALSE, verbose=FALSE)
}

\arguments{
 \item{filename, path}{The filename and path of the generated
       \code{\link{FastaReferenceFile}}.}
 \item{...}{Not used.}
 \item{overwrite}{If \code{\link[base:logical]{TRUE}} and the output file already exists, then
   it is overwritten, otherwise an exception is thrown.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link{FastaReferenceFile}}.
}

\author{Henrik Bengtsson}



\keyword{internal}
\keyword{methods}
\keyword{IO}
