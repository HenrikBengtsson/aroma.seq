%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractAlignment.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{AbstractAlignment}
\docType{class}
\alias{AbstractAlignment}


 \title{The AbstractAlignment class}

 \description{
  Package:  aroma.seq \cr
\bold{Class AbstractAlignment}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\emph{\code{AbstractAlignment}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.seq]{Bowtie2Alignment}, \link[aroma.seq]{BwaAlignment}\cr

public abstract static class \bold{AbstractAlignment}\cr
extends \link[R.oo]{Object}\cr



  ...
 }

 \usage{AbstractAlignment(dataSet=NULL, indexSet=NULL, tags="*", rgSet=NULL, ...)}

 \arguments{
  \item{dataSet}{An \code{\link{FastqDataSet}}.}
  \item{indexSet}{An \code{\link{AbstractIndexSet}}.}
  \item{tags}{Additional tags for the output data sets.}
  \item{rgSet}{(optional) An \code{\link{SamReadGroup}} for added
    SAM read group to the results.}
  \item{...}{Additional alignment arguments.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getIndexSet} \tab  -\cr
 \tab \code{getInputDataSet} \tab  -\cr
 \tab \code{getName} \tab  -\cr
 \tab \code{getOptionalArguments} \tab  -\cr
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{getParameters} \tab  -\cr
 \tab \code{getParametersAsString} \tab  -\cr
 \tab \code{getPath} \tab  -\cr
 \tab \code{getRootPath} \tab  -\cr
 \tab \code{getTags} \tab  -\cr
 \tab \code{nbrOfFiles} \tab  -\cr
 \tab \code{process} \tab  -\cr
 \tab \code{setTags} \tab  -\cr
}


 \bold{Methods inherited from Object}:\cr
asThis, getChecksum, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save


 }

 \author{Henrik Bengtsson}
\keyword{classes}