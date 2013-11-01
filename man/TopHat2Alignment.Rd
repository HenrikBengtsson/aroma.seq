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
\code{~~~~~~~+--}\code{\link[aroma.seq]{AbstractAlignment}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{TopHat2Alignment}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{TopHat2Alignment}\cr
extends \emph{\link[aroma.seq]{AbstractAlignment}}\cr



 ...
}

\usage{
TopHat2Alignment(..., indexSet=NULL, geneModelFile=NULL)
}

\arguments{
 \item{...}{Arguments passed to \code{\link{AbstractAlignment}}.}
 \item{indexSet}{An \code{\link{Bowtie2IndexSet}}.}
 \item{geneModelFile}{Gene model (transcriptome) gtf/gff file.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getOutputDataSet} \tab  -\cr
 \tab \code{process} \tab  -\cr
}


 \bold{Methods inherited from AbstractAlignment}:\cr
as.character, getAsteriskTags, getFlavor, getFullName, getIndexSet, getInputDataSet, getName, getOptionalArguments, getOutputDataSet, getParameterSets, getPath, getRootPath, getTags, isPaired, nbrOfFiles, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\section{Supported operating systems}{
  This method is available on Linux, OSX, and Windows [1].
}

\author{Taku Tokuyasu}

\references{
 [1] TopHat, University of Maryland, 2013.
     \url{http://http://tophat.cbcb.umd.edu/}
}
\keyword{classes}
