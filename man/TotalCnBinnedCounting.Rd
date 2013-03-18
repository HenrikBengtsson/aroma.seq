%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  TotalCnBinnedCounting.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{TotalCnBinnedCounting}
\docType{class}
\alias{TotalCnBinnedCounting}


 \title{The TotalCnBinnedCounting class}

 \description{
  Package:  aroma.seq \cr
\bold{Class TotalCnBinnedCounting}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.cn]{TotalCnSmoothing}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{TotalCnBinnedCounting}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{TotalCnBinnedCounting}\cr
extends \emph{\link[aroma.cn]{TotalCnSmoothing}}\cr



 }

 \usage{
TotalCnBinnedCounting(..., .reqSetClass="BamDataSet")
}

 \arguments{
  \item{...}{Arguments passed to \code{\link[aroma.cn]{TotalCnSmoothing}}.}
  \item{.reqSetClass}{(internal) ...}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\emph{No methods defined}.


 \bold{Methods inherited from TotalCnSmoothing}:\cr
getAsteriskTags, getOutputDataSet0, getOutputFileClass, getOutputFileExtension, getOutputFiles, getOutputFileSetClass, getParameters, getPath, getRootPath, getTargetPositions, getTargetUgpFile, process, smoothRawCopyNumbers

 \bold{Methods inherited from AromaTransform}:\cr
as.character, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameters, getParametersAsString, getParameterSets

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, findRspReportTemplate, asThis


 }

 \author{Henrik Bengtsson}
\keyword{classes}
