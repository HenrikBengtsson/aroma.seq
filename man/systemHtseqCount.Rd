%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  systemHTSeqCount.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{systemHTSeqCount}
\alias{systemHTSeqCount.default}
\alias{systemHTSeqCount}

\title{Wrapper for htseq-count}

\description{
 Wrapper for htseq-count.
}

\usage{
\method{systemHTSeqCount}{default}(commandName="htseq-count", ...,
  system2ArgsList=list(stdout = TRUE, stderr = FALSE), .fake=FALSE, verbose=FALSE)
}

\arguments{
  \item{commandName}{A \code{\link[base]{character}} string specifying the htseq-count command.}
  \item{...}{Additional arguments specifying htseq-count command line switches.}
  \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
  \item{.fake}{If \code{\link[base:logical]{TRUE}}, the executable is not called.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\author{Taku Tokuyasu}


\keyword{internal}
