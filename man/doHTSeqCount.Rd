%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  doHTSeqCount.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{doHTSeqCount}
\alias{doHTSeqCount.default}
\alias{doHTSeqCount}
\alias{doHTSeqCount.BamDataSet}

\title{Counting reads in features}

\description{
 Counting reads in features based on [1].
}

\usage{
  doHTSeqCount(...)
  \method{doHTSeqCount}{BamDataSet}(dataSet, transcripts, ..., verbose=FALSE)
}

\arguments{
 \item{dataSet}{A \code{\link{BamDataSet}}.}
 \item{transcripts}{A \code{\link{GtfDataFile}}.}
 \item{...}{Additional arguments passed to \code{\link{HTSeqCounting}}.}
 \item{verbose}{See \code{\link{Verbose}}.}
}

\value{
  Returns a \code{\link{HTSeqCountDataSet}}.
}

\references{
 [1] Simon Anders, \emph{HTSeq: Analysing high-throughput sequencing
     data with Python}, EMBL, Jan 2014.
     \url{http://www-huber.embl.de/users/anders/HTSeq/} \cr
}

\author{Henrik Bengtsson}

\seealso{
 For more details, see \code{\link{HTSeqCounting}}.
}


\keyword{internal}
