%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SamTransform.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SamTransform}
\docType{class}
\alias{SamTransform}


\title{The SamTransform class}

\description{
 Package:  aroma.seq \cr
\bold{Class SamTransform}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.seq]{AromaSeqTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{SamTransform}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.seq]{BamMerger}, \link[aroma.seq]{PicardDuplicateRemoval}\cr

public abstract static class \bold{SamTransform}\cr
extends \emph{\link[aroma.seq]{AromaSeqTransform}}\cr



 A SamTransform is an \code{\link{AromaSeqTransform}} that takes
 \code{\link{BamDataSet}}:s (or \code{\link{SamDataSet}}:s) as input.
}

\usage{
SamTransform(..., .className="BamDataSet")
}

\arguments{
 \item{...}{Arguments passed to \code{\link{AromaSeqTransform}}.}
 \item{.className}{A \code{\link[base]{character}} string specifying what class
  of data sets to accept.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from AromaSeqTransform}:\cr
as.character, findFilesTodo, getAcronym, getAsteriskTags, getFlavor, getFullName, getInputDataSet, getName, getOptionalArguments, getOrganism, getOutputDataSet, getParameters, getPath, getRootPath, getTags, isDone, length, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\author{Henrik Bengtsson}


\keyword{classes}
\keyword{internal}
