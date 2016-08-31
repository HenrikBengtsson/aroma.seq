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
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.seq]{AromaSeqTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{AbstractAlignment}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.seq]{Bowtie2Alignment}, \link[aroma.seq]{BwaAlignment}, \link[aroma.seq]{TopHat2Alignment}\cr

public abstract static class \bold{AbstractAlignment}\cr
extends \emph{\link[aroma.seq]{AromaSeqTransform}}\cr



 ...
}

\usage{
AbstractAlignment(dataSet=NULL, indexSet=NULL, rgSet=NULL, ...)
}

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
\emph{No methods defined}.


 \bold{Methods inherited from AromaSeqTransform}:\cr
as.character, findFilesTodo, getAcronym, getAsteriskTags, getFlavor, getFullName, getInputDataSet, getName, getOptionalArguments, getOrganism, getOutputDataSet, getParameters, getPath, getRootPath, getTags, isDone, length, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, findRspReportTemplate, asThis


}

\author{Henrik Bengtsson}


\keyword{classes}
\keyword{internal}
