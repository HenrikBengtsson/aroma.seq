%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BamDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{validate.BamDataFile}
\alias{validate.BamDataFile}
\alias{BamDataFile.validate}
\alias{validate,BamDataFile-method}
\alias{validate.SamDataFile}

\title{Validates a BAM (or SAM) file}

\description{
  Validates a BAM (or SAM) file.
}

\usage{
\method{validate}{BamDataFile}(this, method=c("picard"), ..., verbose=FALSE)
}

\arguments{
 \item{method}{A \code{\link[base]{character}} string specifying which validation method to use.}
 \item{...}{Additional arguments passed to the internal validation method.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns (invisibly) what the internal validation method returns,
  or throws an exception if an error is detected.
}

\seealso{
  Internally \code{picardValidateSamFile()} is used.
}

\author{Henrik Bengtsson}


\keyword{internal}
\keyword{methods}
