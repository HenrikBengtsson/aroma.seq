%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BamDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{BamDataFile}
\docType{class}
\alias{BamDataFile}


\title{The abstract BamDataFile class}

\description{
 Package:  aroma.seq \cr
\bold{Class BamDataFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{BamDataFile}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BamDataFile}\cr
extends \emph{\link[R.filesets]{GenericDataFile}}\cr



 A BamDataFile object represents a BAM file.
}

\usage{
BamDataFile(...)
}

\arguments{
  \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFile}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{buildIndex} \tab  -\cr
 \tab \code{convertToSam} \tab  -\cr
 \tab \code{directoryStructure} \tab  -\cr
 \tab \code{doQDNAseq} \tab  -\cr
 \tab \code{extractReadStartPositions} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getHeader} \tab  -\cr
 \tab \code{getIndexFile} \tab  -\cr
 \tab \code{getIndexStats} \tab  -\cr
 \tab \code{getPlatform} \tab  -\cr
 \tab \code{getReadCounts} \tab  -\cr
 \tab \code{getReadGroup} \tab  -\cr
 \tab \code{getReadGroups} \tab  -\cr
 \tab \code{getTargetLengths} \tab  -\cr
 \tab \code{getTargetNames} \tab  -\cr
 \tab \code{getTargets} \tab  -\cr
 \tab \code{getTotalTargetLength} \tab  -\cr
 \tab \code{hasIndex} \tab  -\cr
 \tab \code{isSorted} \tab  -\cr
 \tab \code{nbrOfMappedReads} \tab  -\cr
 \tab \code{nbrOfReads} \tab  -\cr
 \tab \code{nbrOfSeqs} \tab  -\cr
 \tab \code{nbrOfTargets} \tab  -\cr
 \tab \code{nbrOfUnmappedReads} \tab  -\cr
 \tab \code{readDataFrame} \tab  -\cr
 \tab \code{readReadPositions} \tab  -\cr
 \tab \code{replaceAllReadGroups} \tab  -\cr
 \tab \code{sort} \tab  -\cr
 \tab \code{sortByName} \tab  -\cr
 \tab \code{sortByPosition} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFile}:\cr
as.character, clone, compareChecksum, copyTo, equals, fromFile, getAttribute, getAttributes, getChecksum, getChecksumFile, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFileSize, getFileType, getFilename, getFilenameExtension, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, is.na, isFile, isGzipped, linkTo, readChecksum, renameTo, renameToUpperCaseExt, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validate, validateChecksum, writeChecksum, directoryItem, directoryItems, directoryStructure, directoryStructure<-, getDataSetName, getOrganismName, getSampleName, loadCacheFile, loadCacheFileItem, memoizedCall2, saveCacheFile, saveCacheFileItem, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\author{Henrik Bengtsson}

\references{
 [1] The SAM Format Specification Working Group,
     \emph{The SAM Format Specification}, Sept 7, 2011.\cr
}

\seealso{
  An object of this class is typically part of an
  \code{\link{BamDataSet}}.
}
\keyword{classes}
