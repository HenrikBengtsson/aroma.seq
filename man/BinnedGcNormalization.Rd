%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BinnedGcNormalization.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{BinnedGcNormalization}
\docType{class}
\alias{BinnedGcNormalization}


\title{The abstract BinnedGcNormalization class}

\description{
 Package:  aroma.seq \cr
\bold{Class BinnedGcNormalization}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{BinnedGcNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BinnedGcNormalization}\cr
extends \emph{\link[aroma.core]{AromaTransform}}\cr



}

\usage{
BinnedGcNormalization(dataSet=NULL, ..., .reqSetClass="AromaUnitTotalCnBinarySet")
}

\arguments{
 \item{dataSet}{An \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}}.}
 \item{...}{Arguments passed to \code{\link[aroma.core]{AromaTransform}}.}
 \item{.reqSetClass}{(internal only)}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{getAsteriskTags} \tab  -\cr
 \tab \code{getGcContentFile} \tab  -\cr
 \tab \code{getOutputDataSet0} \tab  -\cr
 \tab \code{getOutputFileClass} \tab  -\cr
 \tab \code{getOutputFileExtension} \tab  -\cr
 \tab \code{getOutputFileSetClass} \tab  -\cr
 \tab \code{getOutputFiles} \tab  -\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from AromaTransform}:\cr
as.character, getAsteriskTags, getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getRootPath, getTags, isDone, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, findRspReportTemplate, asThis

 "public"
}

\author{Henrik Bengtsson}
\keyword{classes}
