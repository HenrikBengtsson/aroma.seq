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
\code{~~+--}\code{\link[aroma.core]{AromaTransform}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{BinnedGcNormalization}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BinnedGcNormalization}\cr
extends \emph{\link[aroma.core]{AromaTransform}}\cr



 }

 \usage{BinnedGcNormalization(dataSet=NULL, ..., .reqSetClass="AromaUnitTotalCnBinarySet")}

 \arguments{
  \item{dataSet}{An \code{\link[aroma.core]{AromaUnitTotalCnBinarySet}}.}
  \item{...}{Arguments passed to \code{\link[aroma.core]{AromaTransform}}.}
  \item{.reqSetClass}{(internal only)}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{getGcContentFile} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from AromaTransform}:\cr
getExpectedOutputFiles, getExpectedOutputFullnames, getFullName, getInputDataSet, getName, getOutputDataSet, getOutputDataSet0, getOutputFiles, getPath, getTags, isDone, process, setTags

 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
