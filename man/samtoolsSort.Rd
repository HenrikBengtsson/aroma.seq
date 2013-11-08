%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  samtoolsSort.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{samtoolsSort}
\alias{samtoolsSort.default}
\alias{samtoolsSort}

\title{Calls the samtools 'sort' command}

\description{
 Calls the samtools 'sort' command.
}

\usage{
\method{samtoolsSort}{default}(pathname, pathnameD, ..., verbose=FALSE)
}

\arguments{
  \item{pathname}{A SAM/BAM file.}
  \item{pathnameD}{The destination pathname.}
  \item{...}{Additional arguments specifying samtools 'sort' switches
    passed to \code{\link{systemSamtools}}().}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Henrik Bengtsson}


\keyword{internal}
