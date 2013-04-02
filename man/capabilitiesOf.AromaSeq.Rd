%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.AromaSeq.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AromaSeq$capabilitiesOf}
\alias{AromaSeq$capabilitiesOf}
\alias{capabilitiesOf.AromaSeq}
\alias{AromaSeq.capabilitiesOf}
\alias{capabilitiesOf,AromaSeq-method}
\alias{AromaSeq.isCapableOf}
\alias{isCapableOf.AromaSeq}
\alias{isCapableOf,AromaSeq-method}


\title{Checks which tools are supported}

\description{
  Checks which tools are supported.
}

\usage{
## Static method (use this):
## AromaSeq$capabilitiesOf(what=NULL, ...)

## Don't use the below:
\method{capabilitiesOf}{AromaSeq}(static, what=NULL, ...)
}

\arguments{
 \item{what}{Optional \code{\link[base]{character}} \code{\link[base]{vector}} of which tools to check.}
 \item{...}{Not used.}
}

\value{
  Returns a \code{\link[base]{logical}} named \code{\link[base]{character}} \code{\link[base]{vector}}.
}

\examples{
  # Display which sequencing tools are supported by the package
  print(capabilitiesOf(aroma.seq))

  # Check whether BWA is supported
  print(isCapableOf(aroma.seq, "bwa"))
}

\author{Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
