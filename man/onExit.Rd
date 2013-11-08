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


\keyword{internal}
