%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SamDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{SamDataSet}
\docType{class}
\alias{SamDataSet}


 \title{The SamDataSet class}

 \description{
  Package:  aroma.seq \cr
\bold{Class SamDataSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{SamDataSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{SamDataSet}\cr
extends \link[R.filesets]{GenericDataFileSet}\cr



  An SamDataSet object represents a set of \code{\link{SamDataFile}}:s.
 }

 \usage{SamDataSet(files=NULL, ...)}

 \arguments{
   \item{files}{A \code{\link[base]{list}} of \code{\link{SamDataFile}}:s.}
   \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFileSet}}.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byPath} \tab  -\cr
 \tab \code{\link[aroma.seq:convertToBamDataSet.SamDataSet]{convertToBamDataSet}} \tab Converts the SAM set into a BAM set.\cr
}


 \bold{Methods inherited from GenericDataFileSet}:\cr
getFullNameTranslatorSet, getParentName, append, appendFiles, appendFullNamesTranslator, appendFullNamesTranslatorBydata.frame, appendFullNamesTranslatorByfunction, appendFullNamesTranslatorBylist, appendFullNamesTranslatorByNULL, appendFullNamesTranslatorByTabularTextFile, appendFullNamesTranslatorByTabularTextFileSet, as.list, byName, byPath, clearCache, clearFullNamesTranslator, copyTo, equals, extract, findByName, getChecksum, getChecksumObjects, getDefaultFullName, getFile, getFileClass, getFileSize, getFullNames, getNames, getPath, getPathnames, getSubdirs, hasFile, indexOf, lapply, nbrOfFiles, sapply, seq, setFullNamesTranslator, sortBy, update2, updateFullName, updateFullNames, validate

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, registerFinalizer, save, print


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{classes}
