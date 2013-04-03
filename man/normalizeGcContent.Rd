%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  normalizeGCContent.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{normalizeGcContent}
\alias{normalizeGcContent.default}
\alias{normalizeGcContent}

\title{Normalizes signals for GC-content effects}

\description{
 Normalizes signals for GC-content effects. Some or all signals are used to estimated the
 normalization function.  All signals are normalized.
}

\usage{
\method{normalizeGcContent}{default}(y, gcContent, targetFcn=NULL, subsetToFit=NULL, onMissing=c("ignore", "median"), .isLogged=TRUE, ..., .returnFit=FALSE)
}

\arguments{
  \item{y}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of length K of signals to be normalized.}
  \item{gcContent}{A \code{\link[base]{numeric}} \code{\link[base]{vector}} of GC fractions.}
  \item{targetFcn}{An optional \code{\link[base]{function}}.
    If \code{\link[base]{NULL}}, the data is normalized to have constant GC content effects
    (all equal to zero on the log-scale).}
  \item{subsetToFit}{The subset of data points used to fit the
     normalization function.
     If \code{\link[base]{NULL}}, all data points are considered.}
  \item{onMissing}{Specifies how data points for which there is no
     fragment length is normalized.
     If \code{"ignore"}, the values are not modified.
     If \code{"median"}, the values are updated to have the same
     robust average as the other data points.
  }
  \item{.isLogged}{A \code{\link[base]{logical}}.}
  \item{...}{Additional arguments passed to \code{\link[stats]{lowess}}.}
  \item{.returnFit}{A \code{\link[base]{logical}}.}
}

\value{
  Returns a \code{\link[base]{numeric}} \code{\link[base]{vector}} of the normalized signals.
}

\author{Henrik Bengtsson}




\keyword{nonparametric}
\keyword{robust}
\keyword{internal}
