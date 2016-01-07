%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  AbstractIndexSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{AbstractIndexSet}
\docType{class}
\alias{AbstractIndexSet}

\alias{AbstractIndexSet.getDefaultFilePatterns}
\alias{getDefaultFilePatterns.AbstractIndexSet}
\alias{getDefaultFilePatterns,AbstractIndexSet-method}


\title{The AbstractIndexSet class}

\description{
 Package:  aroma.seq \cr
\bold{Class AbstractIndexSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{AbstractIndexSet}}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.seq]{Bowtie2IndexSet}, \link[aroma.seq]{BwaIndexSet}\cr

public abstract static class \bold{AbstractIndexSet}\cr
extends \link[R.filesets]{GenericDataFileSet}\cr



 An AbstractIndexSet object represents a set of \code{\link{AbstractIndexFile}}:s.
}

\usage{
AbstractIndexSet(files=NULL, ...)
}

\arguments{
  \item{files}{A \code{\link[base]{list}} of \code{\link{AbstractIndexFile}}:s.}
  \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFileSet}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byPrefix} \tab  -\cr
 \tab \code{getFastaReferenceFile} \tab  -\cr
 \tab \code{getIndexPrefix} \tab  -\cr
 \tab \code{getOrganism} \tab  -\cr
 \tab \code{isComplete} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFileSet}:\cr
[, [[, anyDuplicated, anyNA, append, appendFiles, appendFullNamesTranslator, appendFullNamesTranslatorByNULL, appendFullNamesTranslatorByTabularTextFile, appendFullNamesTranslatorByTabularTextFileSet, appendFullNamesTranslatorBydata.frame, appendFullNamesTranslatorByfunction, appendFullNamesTranslatorBylist, as.character, as.list, byName, byPath, c, clearCache, clearFullNamesTranslator, clone, copyTo, dsApply, dsApplyInPairs, duplicated, equals, extract, findByName, findDuplicated, getChecksum, getChecksumFileSet, getChecksumObjects, getDefaultFullName, getFile, getFileClass, getFileSize, getFiles, getFullNames, getNames, getOneFile, getPath, getPathnames, getSubdirs, gunzip, gzip, hasFile, indexOf, is.na, names, nbrOfFiles, rep, resetFullNames, setFullNamesTranslator, sortBy, unique, update2, updateFullName, updateFullNames, validate, directoryItem, directoryItems, directoryStructure, directoryStructure<-, updateDirectoryStructure, getFullNameTranslatorSet, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\author{Henrik Bengtsson}


\keyword{classes}
\keyword{internal}
