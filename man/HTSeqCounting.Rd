%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  HTSeqCounting.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{HTSeqCounting}
\docType{class}
\alias{HTSeqCounting}


\title{The HTSeqCounting class}

\description{
 Package:  aroma.seq \cr
\bold{Class HTSeqCounting}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.seq]{AromaSeqTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{HTSeqCounting}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{HTSeqCounting}\cr
extends \emph{\link[aroma.seq]{AromaSeqTransform}}\cr



 ...
}

\usage{
HTSeqCounting(dataSet=NULL, transcripts=NULL, ...)
}

\arguments{
 \item{dataSet}{A \code{\link{BamDataSet}}.}
 \item{transcripts}{A \code{\link{GtfDataFile}} specifying a gene model.}
 \item{...}{Arguments passed to \code{\link{AbstractAlignment}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from AromaSeqTransform}:\cr
as.character, findFilesTodo, getAcronym, getAsteriskTags, getFlavor, getFullName, getInputDataSet, getName, getOptionalArguments, getOrganism, getOutputDataSet, getParameters, getPath, getRootPath, getTags, isDone, length, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\section{Supported operating systems}{
  ...
}

\author{Henrik Bengtsson}

\references{
 [1] Simon Anders, \emph{HTSeq: Analysing high-throughput sequencing
     data with Python}, EMBL, Jan 2014.
     \url{http://www-huber.embl.de/users/anders/HTSeq/} \cr
}
\keyword{classes}
