%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Bowtie2Alignment.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{Bowtie2Alignment}
\docType{class}
\alias{Bowtie2Alignment}


 \title{The Bowtie2Alignment class}

 \description{
  Package:  aroma.seq \cr
\bold{Class Bowtie2Alignment}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.seq]{AbstractAlignment}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{Bowtie2Alignment}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{Bowtie2Alignment}\cr
extends \emph{\link[aroma.seq]{AbstractAlignment}}\cr



  ...
 }

 \usage{
Bowtie2Alignment(..., indexSet=NULL)
}

 \arguments{
  \item{...}{Arguments passed to \code{\link{AbstractAlignment}}.}
  \item{indexSet}{An \code{\link{Bowtie2IndexSet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getParameters} \tab  -\cr
 \tab \code{\link[aroma.seq:process.Bowtie2Alignment]{process}} \tab Runs the aligner.\cr
}


 \bold{Methods inherited from AbstractAlignment}:\cr
as.character, getAsteriskTags, getFlavor, getFullName, getIndexSet, getInputDataSet, getName, getOptionalArguments, getOutputDataSet, getParameterSets, getPath, getRootPath, getTags, nbrOfFiles, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameters, getParametersAsString, getParameterSets

 \bold{Methods inherited from Object}:\cr
[[, [[<-, $, $<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, findRspReportTemplate, asThis

 "public"
 }

 \author{Henrik Bengtsson and Pierre Neuvial}
\keyword{classes}
