%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastQCReporter.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FastQCReporter}
\docType{class}
\alias{FastQCReporter}


\title{The FastQCReporter class}

\description{
 Package:  aroma.seq \cr
\bold{Class FastQCReporter}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.seq]{AromaSeqTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.seq]{FileGroupsInterface}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{FastQCReporter}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{FastQCReporter}\cr
extends \link[aroma.seq]{FileGroupsInterface}\cr



 A FastQCReporter takes a \code{\link{FastqDataSet}} as input, possibly groups
 the samples, and generates FastQC [1] reports.
 How the grouping is done, can be controlled by a parameter.
}

\usage{
FastQCReporter(dataSet=NULL, groupBy=NULL, ..., .className="FastqDataSet")
}

\arguments{
 \item{dataSet}{An \code{\link{FastqDataSet}}.}
 \item{groupBy}{A \code{\link[base]{character}} string or an explicit named \code{\link[base]{list}},
  specifying which input files should be processed together.}
 \item{...}{Additional arguments passed to \code{\link{AromaSeqTransform}}.}
 \item{.className}{A \code{\link[base]{character}} string specifying what class of
  data sets to accept.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from FileGroupsInterface}:\cr
getGroupBy, getGroupNames, getGroups, getInputDataSet, nbrOfGroups, validateGroups

 \bold{Methods inherited from AromaSeqTransform}:\cr
as.character, findFilesTodo, getAcronym, getAsteriskTags, getFlavor, getFullName, getInputDataSet, getName, getOptionalArguments, getOrganism, getOutputDataSet, getParameters, getPath, getRootPath, getTags, isDone, length, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\section{Supported operating systems}{
  This method is available on Linux, OSX, and Windows [1].
}

\author{Henrik Bengtsson}

\references{
 [1] Simon Andrews,
     FastQC - A quality control tool for high throughput sequence data,
     March 2014.
     \url{http://www.bioinformatics.babraham.ac.uk/projects/fastqc/}
}

\seealso{
  Internally \code{\link[Rsamtools]{mergeBam}} is used for merging, sorting,
  and indexing.
}
\keyword{classes}
