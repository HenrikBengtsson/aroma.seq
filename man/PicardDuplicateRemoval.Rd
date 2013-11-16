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



 This method \emph{flags} reads that are aligned to more than one locus,
 which is done using Picard's 'MarkDuplicates' method [1].

 Note that it is assumed that the input BAM files are already sorted,
 which also means that it can be assumed that the output BAM files
 are sorted.  As with all other methods that outputs BAM files,
 this methods index all outputted BAM files.
}

\usage{
PicardDuplicateRemoval(dataSet=NULL, tags="*", ASSUME_SORTED=TRUE,
  VALIDATION_STRINGENCY="SILENT", ...)
}

\arguments{
 \item{dataSet}{An \code{\link{BamDataSet}}.}
 \item{tags}{Additional tags for the output data sets.}
 \item{ASSUME_SORTED, VALIDATION_STRINGENCY, ...}{
   Additional arguments passed to Picard MarkDuplicates.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{findFilesTodo} \tab  -\cr
 \tab \code{getFullName} \tab  -\cr
 \tab \code{getInputDataSet} \tab  -\cr
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{nbrOfFiles} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\section{Benchmarking}{
 As a very rough guideline, a 1.0GB BAM file takes
 about 10-15 minutes to process using this method.
}

\references{
 [1] Picard, \url{http://picard.sourceforge.net/}\cr
}

\author{Henrik Bengtsson}
\keyword{classes}
