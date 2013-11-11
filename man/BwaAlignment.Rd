%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BwaAlignment.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{BwaAlignment}
\docType{class}
\alias{BwaAlignment}


\title{The BwaAlignment class}

\description{
 Package:  aroma.seq \cr
\bold{Class BwaAlignment}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[aroma.core]{ParametersInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[aroma.seq]{AbstractAlignment}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{BwaAlignment}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BwaAlignment}\cr
extends \emph{\link[aroma.seq]{AbstractAlignment}}\cr



 ...
}

\usage{
BwaAlignment(..., indexSet=NULL)
}

\arguments{
 \item{...}{Arguments passed to \code{\link{AbstractAlignment}}.}
 \item{indexSet}{An \code{\link{BwaIndexSet}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.seq:process.BwaAlignment]{process}} \tab Runs the BWA aligner.\cr
}


 \bold{Methods inherited from AbstractAlignment}:\cr
as.character, findFilesTodo, getAsteriskTags, getFlavor, getFullName, getIndexSet, getInputDataSet, getName, getOptionalArguments, getOrganism, getOutputDataSet, getParameterSets, getPath, getRootPath, getTags, isPaired, nbrOfFiles, process, setTags

 \bold{Methods inherited from ParametersInterface}:\cr
getParameterSets, getParameters, getParametersAsString

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\references{
  [1] Li H. and Durbin R., \emph{Fast and accurate short read alignment
      with Burrows-Wheeler Transform}. Bioinformatics, 2009.\cr
  [2] Li H. and Durbin R., \emph{Fast and accurate long-read alignment
      with Burrows-Wheeler Transform}. Bioinformatics, 2010.\cr
}

\author{Henrik Bengtsson}
\keyword{classes}
