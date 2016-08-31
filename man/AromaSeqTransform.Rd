%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AromaSeqTransform.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AromaSeqTransform}
\docType{class}
\alias{AromaSeqTransform}


\title{The AromaSeqTransform class}

\description{
 Package:  aroma.seq \cr
\bold{Class AromaSeqTransform}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\emph{\code{AromaSeqTransform}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.seq]{AbstractAlignment}}, \link[aroma.seq]{BamDownsampler}, \link[aroma.seq]{BamMerger}, \link[aroma.seq]{Bowtie2Alignment}, \link[aroma.seq]{BwaAlignment}, \link[aroma.seq]{FastQCReporter}, \link[aroma.seq]{FastqDownsampler}, \link[aroma.seq]{HTSeqCounting}, \link[aroma.seq]{PicardDuplicateRemoval}, \link[aroma.seq]{QDNAseqEstimation}, \emph{\link[aroma.seq]{SamTransform}}, \link[aroma.seq]{TopHat2Alignment}\cr

public abstract static class \bold{AromaSeqTransform}\cr
extends \link[aroma.core]{ParametersInterface}\cr



 ...
}

\usage{
AromaSeqTransform(dataSet=NULL, tags="*", flavor=NULL, ...,
  .className="GenericDataFileSet")
}

\arguments{
 \item{dataSet}{An \code{\link{GenericDataFileSet}}.}
 \item{tags}{Tags appended to the output data sets.}
 \item{flavor}{An optional \code{\link[base]{character}} string.}
 \item{...}{Additional named parameters to be recorded.}
 \item{.className}{A \code{\link[base]{character}} string specifying what class of
  data sets to accept.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{length} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, findRspReportTemplate, asThis


}

\author{Henrik Bengtsson}


\keyword{classes}
\keyword{internal}
