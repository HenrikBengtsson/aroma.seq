%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  BamDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{BamDataSet}
\docType{class}
\alias{BamDataSet}


\title{The BamDataSet class}

\description{
 Package:  aroma.seq \cr
\bold{Class BamDataSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.seq]{AromaPathnameInterface}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{\link[aroma.seq]{AromaSeqDataFileSet}}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~~~~~~+--}\code{BamDataSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{BamDataSet}\cr
extends \link[aroma.seq]{AromaSeqDataFileSet}\cr



 An BamDataSet object represents a set of \code{\link{BamDataFile}}:s.
}

\usage{
BamDataSet(files=NULL, ...)
}

\arguments{
  \item{files}{A \code{\link[base]{list}} of \code{\link{BamDataFile}}:s.}
  \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFileSet}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{byName} \tab  -\cr
 \tab \code{byPath} \tab  -\cr
 \tab \code{convertToSam} \tab  -\cr
 \tab \code{countNucleotides} \tab  -\cr
 \tab \code{doDownsample} \tab  -\cr
 \tab \code{doHTSeqCount} \tab  -\cr
 \tab \code{doQDNAseq} \tab  -\cr
 \tab \code{findByName} \tab  -\cr
 \tab \code{getChipType} \tab  -\cr
 \tab \code{getPlatform} \tab  -\cr
 \tab \code{isPaired} \tab  -\cr
 \tab \code{mpileup} \tab  -\cr
 \tab \code{sortByName} \tab  -\cr
 \tab \code{sortByPosition} \tab  -\cr
 \tab \code{splitByTargetType} \tab  -\cr
}


 \bold{Methods inherited from AromaSeqDataFileSet}:\cr
doBWA, getDefaultFullName, isCompatibleWith

 \bold{Methods inherited from AromaPathnameInterface}:\cr
directoryStructure, getOrganism

 \bold{Methods inherited from GenericDataFileSet}:\cr
[, [[, anyDuplicated, anyNA, append, appendFiles, appendFullNamesTranslator, appendFullNamesTranslatorByNULL, appendFullNamesTranslatorByTabularTextFile, appendFullNamesTranslatorByTabularTextFileSet, appendFullNamesTranslatorBydata.frame, appendFullNamesTranslatorByfunction, appendFullNamesTranslatorBylist, as.character, as.list, byName, byPath, c, clearCache, clearFullNamesTranslator, clone, copyTo, dsApplyInPairs, duplicated, equals, extract, findByName, findDuplicated, getChecksum, getChecksumFileSet, getChecksumObjects, getDefaultFullName, getFile, getFileClass, getFileSize, getFiles, getFullNames, getNames, getOneFile, getPath, getPathnames, getSubdirs, gunzip, gzip, hasFile, indexOf, is.na, names, nbrOfFiles, rep, resetFullNames, setFullNamesTranslator, sortBy, unique, update2, updateFullName, updateFullNames, validate, directoryItem, directoryItems, directoryStructure, directoryStructure<-, updateDirectoryStructure, getFullNameTranslatorSet, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, names, objectSize, print, save, findRspReportTemplate, asThis


}

\author{Henrik Bengtsson}
\keyword{classes}
