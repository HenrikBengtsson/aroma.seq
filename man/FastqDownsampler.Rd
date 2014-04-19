%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastqDownsampler.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FastqDownsampler}
\docType{class}
\alias{FastqDownsampler}


\title{The FastqDownsampler class}

\description{
 Package:  aroma.seq \cr
\bold{Class FastqDownsampler}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.seq]{AromaSeqTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{FastqDownsampler}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{FastqDownsampler}\cr
extends \emph{\link[aroma.seq]{AromaSeqTransform}}\cr



 ...
}

\usage{
FastqDownsampler(dataSet=NULL, subset=1e+06, seed=NULL, ...)
}

\arguments{
 \item{dataSet}{An \code{\link{FastqDataSet}}.}
 \item{subset}{An \code{\link[base]{integer}} specifying the total number of reads to sample,
   or a \code{\link[base]{double}} specifying the fraction of total number of reads to sample.}
 \item{seed}{An (optional) \code{\link[base]{integer}} specifying the random seed to be
    set before sampling indices.  The random seed is set to its original
    state when exiting.  If \code{\link[base]{NULL}}, it is not set.}
 \item{...}{Additional arguments passed to \code{\link{AromaSeqTransform}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from AromaSeqTransform}:\cr
as.character, findFilesTodo, getAcronym, getAsteriskTags, getFlavor, getFullName, getInputDataSet, getName, getOptionalArguments, getOrganism, getOutputDataSet, getParameters, getPath, getRootPath, getTags, isDone, length, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\seealso{
 Internally, the \code{\link[ShortRead]{FastqSampler}} method is used.
}

\author{Henrik Bengtsson}
\keyword{classes}
