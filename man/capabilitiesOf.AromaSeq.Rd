%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  999.AromaSeq.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{capabilitiesOf.AromaSeq}
\alias{capabilitiesOf.AromaSeq}
\alias{AromaSeq.capabilitiesOf}
\alias{capabilitiesOf.AromaSeq}
\alias{capabilitiesOf,AromaSeq-method}

 \alias{AromaSeq.isCapableOf}
\alias{isCapableOf.AromaSeq}
\alias{isCapableOf,AromaSeq-method}


 \title{Checks which tools are supported}

 \description{
   Checks which tools are supported.
 }

 \synopsis{capabilitiesOf.AromaSeq(static, what=NULL, ...)}
\usage{AromaSeq$capabilitiesOf(what=NULL, ...)}

 \arguments{
  \item{what}{Optional \code{\link[base]{character}} \code{\link[base]{vector}} of which tools to check.}
  \item{...}{Not used.}
 }

 \value{
   Returns a \code{\link[base]{logical}} named \code{\link[base]{character}} \code{\link[base]{vector}}.
 }

 \examples{
   # Display which sequencing tools are supported by the package
   print(capabilitiesOf(aroma.seq))

   # Check whether BWA is supported
   print(isCapableOf(aroma.seq, "bwa"))
 }

 \author{TBA}
\keyword{internal}
\keyword{methods}
