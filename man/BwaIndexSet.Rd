%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BwaIndexSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{BwaIndexSet}
\docType{class}
\alias{BwaIndexSet}

\alias{BwaIndexSet.as.character}
\alias{as.character.BwaIndexSet}
\alias{as.character,BwaIndexSet-method}

\alias{BwaIndexSet.getDefaultFilePatterns}
\alias{getDefaultFilePatterns.BwaIndexSet}
\alias{getDefaultFilePatterns,BwaIndexSet-method}

\alias{BwaIndexSet.getSeqLengths}
\alias{getSeqLengths.BwaIndexSet}
\alias{getSeqLengths,BwaIndexSet-method}

\alias{BwaIndexSet.getSeqNames}
\alias{getSeqNames.BwaIndexSet}
\alias{getSeqNames,BwaIndexSet-method}

\alias{BwaIndexSet.isCompatibleWith}
\alias{isCompatibleWith.BwaIndexSet}
\alias{isCompatibleWith,BwaIndexSet-method}

\alias{BwaIndexSet.isComplete}
\alias{isComplete.BwaIndexSet}
\alias{isComplete,BwaIndexSet-method}

\alias{BwaIndexSet.readAnnData}
\alias{readAnnData.BwaIndexSet}
\alias{readAnnData,BwaIndexSet-method}


\title{The BwaIndexSet class}

\description{
 Package:  aroma.seq \cr
\bold{Class BwaIndexSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.seq]{AbstractIndexSet}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.seq]{SequenceContigsInterface}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{BwaIndexSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BwaIndexSet}\cr
extends \emph{\link[aroma.seq]{SequenceContigsInterface}}\cr



 An BwaIndexSet object represents a set of \code{\link{BwaIndexFile}}:s.
}

\usage{
BwaIndexSet(...)
}

\arguments{
  \item{...}{Arguments passed to \code{\link{AbstractIndexSet}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{getSeqLengths} \tab  -\cr
 \tab \code{isCompatibleWith} \tab  -\cr
 \tab \code{isComplete} \tab  -\cr
 \tab \code{readAnnData} \tab  -\cr
}


 \bold{Methods inherited from SequenceContigsInterface}:\cr
cleanSeqNames, getSeqGenericSummary, getSeqLengths, getSeqNames, getSeqOrdering, getTotalSeqLength, hasDuplicatedSeqs, hasSeqLengths, isCompatibleWithBySeqLengths, isCompatibleWithBySeqNames, isCompatibleWithBySeqs, nbrOfSeqs

 \bold{Methods inherited from AbstractIndexSet}:\cr
as.character, byPrefix, getDefaultFilePatterns, getFastaReferenceFile, getIndexPrefix, getOrganism, isComplete

 \bold{Methods inherited from GenericDataFileSet}:\cr
[, [[, anyDuplicated, anyNA, append, appendFiles, appendFullNamesTranslator, appendFullNamesTranslatorByNULL, appendFullNamesTranslatorByTabularTextFile, appendFullNamesTranslatorByTabularTextFileSet, appendFullNamesTranslatorBydata.frame, appendFullNamesTranslatorByfunction, appendFullNamesTranslatorBylist, as.character, as.list, byName, byPath, c, clearCache, clearFullNamesTranslator, clone, copyTo, dsApplyInPairs, duplicated, equals, extract, findByName, findDuplicated, getChecksum, getChecksumFileSet, getChecksumObjects, getDefaultFullName, getFile, getFileClass, getFileSize, getFiles, getFullNames, getNames, getOneFile, getPath, getPathnames, getSubdirs, gunzip, gzip, hasFile, indexOf, is.na, names, nbrOfFiles, rep, resetFullNames, setFullNamesTranslator, sortBy, unique, update2, updateFullName, updateFullNames, validate, directoryItem, directoryItems, directoryStructure, directoryStructure<-, updateDirectoryStructure, getFullNameTranslatorSet, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, findRspReportTemplate, asThis


}

\author{Henrik Bengtsson}


\keyword{classes}
\keyword{internal}
