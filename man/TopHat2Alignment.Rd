%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  TopHat2Alignment.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{TopHat2Alignment}
\docType{class}
\alias{TopHat2Alignment}


\title{The TopHat2Alignment class}

\description{
 Package:  aroma.seq \cr
\bold{Class TopHat2Alignment}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.seq]{AromaSeqTransform}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.seq]{AbstractAlignment}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.seq]{FileGroupsInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{TopHat2Alignment}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{TopHat2Alignment}\cr
extends \link[aroma.seq]{FileGroupsInterface}\cr



 ...
}

\usage{
TopHat2Alignment(..., groupBy=NULL, indexSet=NULL, transcripts=NULL)
}

\arguments{
 \item{...}{Arguments passed to \code{\link{AbstractAlignment}}.}
 \item{groupBy}{A \code{\link[base]{character}} string or an explicit named \code{\link[base]{list}},
  specifying which input files should be processed together.}
 \item{indexSet}{An \code{\link{Bowtie2IndexSet}}.}
 \item{transcripts}{Gene model (transcriptome) GTF/GFF3 file.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from FileGroupsInterface}:\cr
getGroupBy, getGroupNames, getGroups, getInputDataSet, nbrOfGroups, validateGroups

 \bold{Methods inherited from AbstractAlignment}:\cr
as.character, getAcronym, getAsteriskTags, getIndexSet, getOrganism, getOutputDataSet, getRootPath, isPaired

 \bold{Methods inherited from AromaSeqTransform}:\cr
as.character, findFilesTodo, getAcronym, getAsteriskTags, getFlavor, getFullName, getInputDataSet, getName, getOptionalArguments, getOrganism, getOutputDataSet, getParameters, getPath, getRootPath, getTags, isDone, length, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\section{Supported operating systems}{
  This method is available on Linux and OSX [1].
}

\author{Taku Tokuyasu}

\references{
 [1] TopHat, University of Maryland, 2013.
     \url{http://http://tophat.cbcb.umd.edu/} \cr
 [2] Trapnell et al. \emph{Differential gene and transcript expression
     analysis of RNA-seq experiments with TopHat and Cufflinks}.
     Nat Protoc, 2012.\cr
}
\keyword{classes}
