%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  HTSeqCountDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{readDGE.HTSeqCountDataSet}
\alias{readDGE.HTSeqCountDataSet}
\alias{HTSeqCountDataSet.readDGE}
\alias{readDGE,HTSeqCountDataSet-method}
\alias{readDGE}
\alias{HTSeqCountDataSet.extractDGEList}
\alias{extractDGEList.HTSeqCountDataSet}
\alias{extractDGEList,HTSeqCountDataSet-method}

\alias{extractDGEList}

\title{Reads all digital gene expression (DGE) data}

\description{
  Reads all digital gene expression (DGE) data.
}

\usage{
\method{readDGE}{HTSeqCountDataSet}(this, labels=getFullNames(this), ...)
}

\arguments{
 \item{labels}{A \code{\link[base]{character}} string specifying sample names.}
 \item{...}{Additional arguments passed to \code{\link[edgeR]{readDGE}}.}
}

\value{
  Returns a \code{\link[edgeR]{DGEList}} object.
}

\author{Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
