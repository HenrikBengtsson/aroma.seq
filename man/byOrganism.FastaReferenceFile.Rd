%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastaReferenceFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FastaReferenceFile$byOrganism}
\alias{FastaReferenceFile$byOrganism}
\alias{byOrganism.FastaReferenceFile}
\alias{FastaReferenceFile.byOrganism}
\alias{byOrganism,FastaReferenceFile-method}
\alias{FastaReferenceFile.findByOrganism}
\alias{findByOrganism.FastaReferenceFile}
\alias{findByOrganism,FastaReferenceFile-method}


\title{Locates a FASTA file by organism}

\description{
  Locates a FASTA file by organism.
}

\usage{
## Static method (use this):
## FastaReferenceFile$byOrganism(organism, ...)

## Don't use the below:
\method{byOrganism}{FastaReferenceFile}(static, organism, ...)
}

\arguments{
 \item{organism}{A \code{\link[base]{character}} string specifying for which organism a
   file should be retrieved.}
 \item{tags}{(not used) A \code{\link[base]{character}} \code{\link[base]{vector}}.}
 \item{prefix}{(optional) A \code{\link[base]{character}} string specifying an optional
   regular expression prefix to be prepended to \code{pattern} when
   searching for the file.}
 \item{pattern}{A \code{\link[base]{character}} string specifying a regular expression for
   the file to be located.}
 \item{...}{Additional arguments passed to the constructor of
   \code{\link{FastaReferenceFile}} when instantiating the object.}
}

\value{
  Returns a \code{\link{FastaReferenceFile}}.
}

\seealso{
  For more information see \code{\link{FastaReferenceFile}}.
}

\author{Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
