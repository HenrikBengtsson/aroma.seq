%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  PicardDuplicateRemoval.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{PicardDuplicateRemoval}
\docType{class}
\alias{PicardDuplicateRemoval}


 \title{The PicardDuplicateRemoval class}

 \description{
  Package:  aroma.seq \cr
\bold{Class PicardDuplicateRemoval}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{PicardDuplicateRemoval}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{PicardDuplicateRemoval}\cr
extends \link[aroma.core]{ParametersInterface}\cr



  ...
 }

 \usage{PicardDuplicateRemoval(dataSet=NULL, tags="*", ASSUME_SORTED=TRUE, VALIDATION_STRINGENCY="SILENT", ...)}

 \arguments{
  \item{dataSet}{An \code{\link{BamDataSet}}.}
  \item{tags}{Additional tags for the output data sets.}
  \item{ASSUME_SORTED, VALIDATION_STRINGENCY, ...}{
    Additional Picard MarkDuplicates arguments.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getInputDataSet} \tab  -\cr
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{nbrOfFiles} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from ParametersInterface}:\cr
getParameters, getParametersAsString, getParameterSets

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, findRspReportTemplate, asThis


 }

 \author{Henrik Bengtsson}
\keyword{classes}
