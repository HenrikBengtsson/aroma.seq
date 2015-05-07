%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastaReferenceIndexFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FastaReferenceIndexFile}
\docType{class}
\alias{FastaReferenceIndexFile}


\title{The FastaReferenceIndexFile class}

\description{
 Package:  aroma.seq \cr
\bold{Class FastaReferenceIndexFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.seq]{SequenceContigsInterface}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\emph{\code{FastaReferenceIndexFile}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{FastaReferenceIndexFile}\cr
extends \emph{\link[aroma.seq]{SequenceContigsInterface}}\cr



 A FastaReferenceIndexFile object represents a
 FASTA reference index file (FAI).
}

\usage{
FastaReferenceIndexFile(...)
}

\arguments{
  \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFile}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getSeqLengths} \tab  -\cr
}


 \bold{Methods inherited from SequenceContigsInterface}:\cr
cleanSeqNames, getSeqGenericSummary, getSeqLengths, getSeqNames, getSeqOrdering, getTotalSeqLength, hasDuplicatedSeqs, hasSeqLengths, isCompatibleWithBySeqLengths, isCompatibleWithBySeqNames, isCompatibleWithBySeqs, nbrOfSeqs

 \bold{Methods inherited from GenericDataFile}:\cr
as.character, clone, compareChecksum, copyTo, equals, fromFile, getAttribute, getAttributes, getChecksum, getChecksumFile, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFileSize, getFileType, getFilename, getFilenameExtension, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, is.na, isFile, isGzipped, linkTo, readChecksum, renameTo, renameToUpperCaseExt, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validate, validateChecksum, writeChecksum, directoryItem, directoryItems, directoryStructure, directoryStructure<-, getDataSetName, getOrganismName, getSampleName, loadCacheFile, loadCacheFileItem, memoizedCall2, saveCacheFile, saveCacheFileItem, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\author{Henrik Bengtsson}

\seealso{
  [1] \url{http://www.wikipedia.org/wiki/FASTA_format}
}
\keyword{classes}
