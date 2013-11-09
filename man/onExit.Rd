%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  onExit.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{onExit}
\alias{onExit}


\title{Records expressions to be executed when the current function exits}

\description{
 Records expressions to be executed when the current function exits (regardless of cause).
 This methtod differ from \code{\link[base]{on.exit}}() in that it adds the
 option to \emph{prepend} an expression to already existing expressed,
 as an alternative to the default \emph{append}.

 \emph{WARNING: This function is still under development [1,2]. /HB 2013-11-08}
}

\usage{
onExit(expr=NULL, where=c("first", "last", "replace"), ..., exits, envir=parent.frame())
}

\arguments{
  \item{expr}{An \code{\link[base]{expression}} to be recorded.}
  \item{where}{A \code{\link[base]{character}} string specifying whether the \code{\link[base]{expression}}
   should be prepended (\code{"first"}), appended (\code{"last"}), or
   replace already recorded \code{\link[base]{expression}}s (\code{"replace"}).}
  \item{...}{Not used.}
  \item{exits}{(optional) Already recorded.}
  \item{envir}{The \code{\link[base]{environment}} in which \code{\link[base]{expression}}s are recorded.}
}

\author{Henrik Bengtsson}

\references{
 [1] R-devel thread 'WISHLIST: on.exit(..., add=TRUE, where="first")
     to address common use cases', 2013-11-03.
     \url{https://stat.ethz.ch/pipermail/r-devel/2013-November/067876.html}
     \cr
 [2] R-devel thread 'on.exit() & sys.on.exit(): Calling them via eval() does not	work as hoped', 2013-11-03.
     \url{https://stat.ethz.ch/pipermail/r-devel/2013-November/067866.html}
     \cr
}


\keyword{internal}
