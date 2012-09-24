%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.AromaSeq.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{isCapableOf.AromaSeq}
\alias{isCapableOf.AromaSeq}
\alias{AromaSeq.isCapableOf}
\alias{isCapableOf.AromaSeq}
\alias{isCapableOf,AromaSeq-method}


 \title{Checks which tools are supported}

 \description{
   Checks which tools are supported.
 }

 \synopsis{isCapableOf.AromaSeq(static, what=NULL, ...)}
\usage{AromaSeq$isCapableOf(what=NULL, ...)}

 \arguments{
  \item{what}{Optional \code{\link[base]{character}} \code{\link[base]{vector}} of which tools to check.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{logical}} named \code{\link[base]{character}} \code{\link[base]{vector}}.
 }

 \examples{
   # Display which sequencing tools are supported
   print(AromaSeq$isCapableOf())

   # Check whether BWA is supported
   print(AromaSeq$isCapableOf("bwa"))
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
