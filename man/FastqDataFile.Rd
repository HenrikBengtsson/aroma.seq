%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastqDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FastqDataFile}
\docType{class}
\alias{FastqDataFile}


\title{The abstract FastqDataFile class}

\description{
 Package:  aroma.seq \cr
\bold{Class FastqDataFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.seq]{AromaPathnameInterface}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.seq]{AromaSeqDataFile}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{FastqDataFile}}\cr

\bold{Directly known subclasses:}\cr
\emph{\link[aroma.seq]{IlluminaFastqDataFile}}\cr

public abstract static class \bold{FastqDataFile}\cr
extends \link[aroma.seq]{AromaSeqDataFile}\cr



 A FastqDataFile object represents a FASTQ data file.
}

\usage{
FastqDataFile(..., paired=FALSE)
}

\arguments{
  \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFile}}.}
  \item{paired}{If \code{\link[base:logical]{TRUE}}, the data set contains paired-end reads,
    otherwise not.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{getCommonSeqWidth} \tab  -\cr
 \tab \code{getDefaultSamReadGroup} \tab  -\cr
 \tab \code{getGeometry} \tab  -\cr
 \tab \code{getSamReadGroup} \tab  -\cr
 \tab \code{nbrOfSeqs} \tab  -\cr
 \tab \code{readSeqFile} \tab  -\cr
 \tab \code{report} \tab  -\cr
 \tab \code{setSamReadGroup} \tab  -\cr
}


 \bold{Methods inherited from AromaSeqDataFile}:\cr
getDefaultFullName

 \bold{Methods inherited from AromaPathnameInterface}:\cr
directoryStructure, getOrganism

 \bold{Methods inherited from GenericDataFile}:\cr
as.character, clone, compareChecksum, copyTo, equals, fromFile, getAttribute, getAttributes, getChecksum, getChecksumFile, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFileSize, getFileType, getFilename, getFilenameExtension, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, is.na, isFile, isGzipped, linkTo, readChecksum, renameTo, renameToUpperCaseExt, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validate, validateChecksum, writeChecksum, directoryItem, directoryItems, directoryStructure, directoryStructure<-, getDataSetName, getOrganismName, getSampleName, loadCacheFile, loadCacheFileItem, memoizedCall2, saveCacheFile, saveCacheFileItem, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\author{Henrik Bengtsson}

\references{
  Wikipedia, FASTQ format,
  \url{http://en.wikipedia.org/wiki/FASTQ_format}.\cr
}

\seealso{
  An object of this class is typically part of a \code{\link{FastqDataSet}}.
}
\keyword{classes}
