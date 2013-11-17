%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findFilesTodo.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findFilesTodo}
\alias{findFilesTodo}

\alias{findFilesTodo.AromaSeqTransform}
\alias{findFilesTodo.TotalCnBinnedCounting}

\title{Identifies which files are not yet processed}

\description{
  Identifies which files are not yet processed.
}

\usage{
 \method{findFilesTodo}{AromaSeqTransform}(this, force=FALSE, ...)
 \method{findFilesTodo}{TotalCnBinnedCounting}(this, ...)
}

\arguments{
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{vector}} of named indices.
}

\seealso{
  Internally \code{getOutputDataSet(..., onMissing="NA")} is utilized.
}

\author{Henrik Bengtsson}


\keyword{methods}
\keyword{internal}
