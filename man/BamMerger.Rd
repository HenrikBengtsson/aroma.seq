%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BamMerger.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{BamMerger}
\docType{class}
\alias{BamMerger}


\title{The BamMerger class}

\description{
 Package:  aroma.seq \cr
\bold{Class BamMerger}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.seq]{AromaSeqTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.seq]{SamTransform}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{BamMerger}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BamMerger}\cr
extends \emph{\link[aroma.seq]{SamTransform}}\cr



 A BamMerger takes a \code{\link{BamDataSet}} as input and merges subsets of
 \code{\link{BamDataFile}}:s into single ones.  How the grouping is done, can
 be controlled by a parameter.
}

\usage{
BamMerger(dataSet=NULL, groupBy=NULL, ...)
}

\arguments{
 \item{dataSet}{An \code{\link{BamDataSet}}.}
 \item{groupBy}{A \code{\link[base]{character}} string or an explicit named \code{\link[base]{list}},
  specifying which BAM files should be merged together.}
 \item{...}{Not used.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getGroups} \tab  -\cr
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{nbrOfGroups} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from SamTransform}:\cr
getRootPath

 \bold{Methods inherited from AromaSeqTransform}:\cr
as.character, findFilesTodo, getAcronym, getAsteriskTags, getFlavor, getFullName, getInputDataSet, getName, getOptionalArguments, getOrganism, getOutputDataSet, getParameters, getPath, getRootPath, getTags, isDone, length, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\section{Supported operating systems}{
  This method is available on Linux, OSX, and Windows [1].
}

\author{Henrik Bengtsson}

\seealso{
  Internally \code{\link[Rsamtools]{mergeBam}} is used for merging, sorting,
  and indexing.
}
\keyword{classes}
