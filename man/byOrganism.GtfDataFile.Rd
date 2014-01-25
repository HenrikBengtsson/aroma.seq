%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  GtfDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{GtfDataFile$byOrganism}
\alias{GtfDataFile$byOrganism}
\alias{byOrganism.GtfDataFile}
\alias{GtfDataFile.byOrganism}
\alias{byOrganism,GtfDataFile-method}
\alias{GtfDataFile.findByOrganism}
\alias{findByOrganism.GtfDataFile}
\alias{findByOrganism,GtfDataFile-method}


\title{Locates a GTF file by organism}

\description{
  Locates a GTF file by organism.
}

\usage{
## Static method (use this):
## GtfDataFile$byOrganism(organism, ...)

## Don't use the below:
\method{byOrganism}{GtfDataFile}(static, organism, ...)
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
   \code{\link{GtfDataFile}} when instantiating the object.}
}

\value{
  Returns a \code{\link{GtfDataFile}}.
}

\seealso{
  For more information see \code{\link{GtfDataFile}}.
}

\author{Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
