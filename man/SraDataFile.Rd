%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SraDataFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SraDataFile}
\docType{class}
\alias{SraDataFile}


\title{The abstract SraDataFile class}

\description{
 Package:  aroma.seq \cr
\bold{Class SraDataFile}\cr

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
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\emph{\code{SraDataFile}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{SraDataFile}\cr
extends \link[aroma.seq]{AromaSeqDataFile}\cr



 A SraDataFile object represents a Sequence Read Archive (SRA) file [1].
}

\usage{
SraDataFile(...)
}

\arguments{
  \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFile}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{fastqDump} \tab  -\cr
}


 \bold{Methods inherited from AromaSeqDataFile}:\cr
getDefaultFullName, isCompatibleWith

 \bold{Methods inherited from AromaPathnameInterface}:\cr
directoryStructure, getOrganism

 \bold{Methods inherited from GenericDataFile}:\cr
as.character, clone, compareChecksum, copyTo, equals, fromFile, getAttribute, getAttributes, getChecksum, getChecksumFile, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFileSize, getFileType, getFilename, getFilenameExtension, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, is.na, isFile, isGzipped, linkTo, readChecksum, renameTo, renameToUpperCaseExt, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validate, validateChecksum, writeChecksum, directoryItem, directoryItems, directoryStructure, directoryStructure<-, getDataSetName, getOrganismName, getSampleName, loadCacheFile, loadCacheFileItem, memoizedCall2, saveCacheFile, saveCacheFileItem, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\references{
 [1] \emph{The Sequence Read Archive (SRA)},
     \url{http://www.ncbi.nlm.nih.gov/sra/}, 2014.\cr
}

\seealso{
  An object of this class is typically part of an
  \code{\link{SraDataSet}}.
}

\author{Henrik Bengtsson}
\keyword{classes}
