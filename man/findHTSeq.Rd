%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  findHTSeq.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{findHTSeq}
\alias{findHTSeq}


\title{Locates the HTSeq executable}

\description{
 Locates the HTSeq executable on the current system.
}

\usage{
findHTSeq(..., command=c("htseq-count", "htseq-qa"))
}

\arguments{
  \item{...}{Arguments passed to \code{\link{findExternal}}().}
  \item{command}{The type of HTSeq executable to locate.}
}

\details{
 The HTSeq tool is search for as follows:
 \enumerate{
  \item ...
 }
}

\author{Henrik Bengtsson, Taku Tokuyasu}

\seealso{
 [1] HTSeq: Analysing high-throughput sequencing data with Python,
     June 2013.
     \url{http://www-huber.embl.de/users/anders/HTSeq/}
}


\keyword{internal}
