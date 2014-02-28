%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  QDNAseqEstimation.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{QDNAseqEstimation}
\docType{class}
\alias{QDNAseqEstimation}


\title{The QDNAseqEstimation class}

\description{
 Package:  aroma.seq \cr
\bold{Class QDNAseqEstimation}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.seq]{AromaSeqTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{QDNAseqEstimation}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{QDNAseqEstimation}\cr
extends \emph{\link[aroma.seq]{AromaSeqTransform}}\cr


}

\usage{
QDNAseqEstimation(dataSet=NULL, binWidth=NULL, log=TRUE, mappability=50, blacklist=0,
  residual=2, bases=0, ...)
}

\arguments{
 \item{dataSet}{An \code{\link{BamDataSet}}.}
 \item{binWidth}{A positive \code{\link[base]{numeric}} specifying the bin width (in units
   of kbp).
   Alternatively, a \code{\link[Biobase]{AnnotatedDataFrame}} specifying the bins.}
 \item{log}{If \code{\link[base:logical]{TRUE}}, the copy numbers are calculated on the log2 scale.}
 \item{mappability, blacklist, residual, bases}{Post-filter arguments.}
 \item{...}{Additional arguments passed to the QDNAseq method.}
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

\section{Benchmarking}{
  TBA.
}

\references{
 [1] TBA.
}

\author{Henrik Bengtsson}
\keyword{classes}
