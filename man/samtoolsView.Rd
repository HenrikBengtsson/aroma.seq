%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  samtoolsView.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{samtoolsView}
\alias{samtoolsView.default}
\alias{samtoolsView}

\title{Calls the samtools 'view' command}

\description{
 Calls the samtools 'view' command.
}

\usage{
\method{samtoolsView}{default}(pathname, pathnameD, ..., verbose=FALSE)
}

\arguments{
  \item{pathname}{A SAM/BAM file.}
  \item{pathnameD}{The destination pathname.}
  \item{...}{Additional arguments specifying samtools 'view' switches
    passed to \code{\link{systemSamtools}}().}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Henrik Bengtsson}

