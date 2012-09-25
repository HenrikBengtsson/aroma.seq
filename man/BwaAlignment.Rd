%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BwaAlignment.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{BwaAlignment}
\docType{class}
\alias{BwaAlignment}


 \title{The BwaAlignment class}

 \description{
  Package:  aroma.seq \cr
\bold{Class BwaAlignment}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{BwaAlignment}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BwaAlignment}\cr
extends \link[R.oo]{Object}\cr



  ...
 }

 \usage{BwaAlignment(dataSet=NULL, indexSet=NULL, tags="*", ...)}

 \arguments{
  \item{dataSet}{An \code{\link{FastqDataSet}}.}
  \item{indexSet}{An \code{\link{BwaIndexSet}}.}
  \item{tags}{Additional tags for the output data sets.}
  \item{...}{Additional arguments passed to the aligner.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getIndexSet} \tab  -\cr
 \tab \code{getInputDataSet} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{nbrOfFiles} \tab  -\cr
 \tab \code{\link[aroma.seq:process.BwaAlignment]{process}} \tab Runs the BWA aligner.\cr
 \tab \code{setTags} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson and Pierre Neuvial}
\keyword{classes}
