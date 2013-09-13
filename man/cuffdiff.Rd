%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  cuffdiff.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{cuffdiff}
\alias{cuffdiff.default}
\alias{cuffdiff}

\title{Calls cuffdiff on input bam file(s)}

\description{
 Calls cuffdiff on input bam file(s).
}

\usage{
\method{cuffdiff}{default}(transcriptsGtf, bams=NULL, optionsVec=NULL, ..., system2ArgsList=NULL,
  commandName="cuffdiff", verbose=FALSE)
}

\arguments{
  \item{transcriptsGtf}{Gtf file of transcripts ('gene model', e.g. from TopHat)}
  \item{bams}{Vector of pathnames (.sam or .bam)}
  \item{optionsVec}{Vector of named options}
  \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
  \item{...}{...}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
  \item{commandName}{Name of executable}
}

\author{Taku Tokuyasu}


\keyword{internal}
