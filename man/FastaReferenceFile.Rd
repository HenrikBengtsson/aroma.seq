%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastaReferenceFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FastaReferenceFile}
\docType{class}
\alias{FastaReferenceFile}


\title{The FastaReferenceFile class}

\description{
 Package:  aroma.seq \cr
\bold{Class FastaReferenceFile}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFile}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\emph{\code{FastaReferenceFile}}\cr

\bold{Directly known subclasses:}\cr
\cr

public abstract static class \bold{FastaReferenceFile}\cr
extends \emph{\link[R.filesets]{GenericDataFile}}\cr



 A FastaReferenceFile object represents a FASTA reference file.
}

\usage{
FastaReferenceFile(...)
}

\arguments{
  \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFile}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{buildBowtie2IndexSet} \tab  -\cr
 \tab \code{buildBwaIndexSet} \tab  -\cr
 \tab \code{buildDictionary} \tab  -\cr
 \tab \code{buildIndex} \tab  -\cr
 \tab \code{byOrganism} \tab  -\cr
 \tab \code{getIndexFile} \tab  -\cr
 \tab \code{getOrganism} \tab  -\cr
 \tab \code{getSeqLengths} \tab  -\cr
 \tab \code{getSeqNames} \tab  -\cr
 \tab \code{getTotalSeqLengths} \tab  -\cr
 \tab \code{hasIndex} \tab  -\cr
 \tab \code{isCompatibleWith} \tab  -\cr
 \tab \code{nbrOfSeqs} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFile}:\cr
as.character, clone, compareChecksum, copyTo, equals, fromFile, getAttribute, getAttributes, getChecksum, getChecksumFile, getCreatedOn, getDefaultFullName, getExtension, getExtensionPattern, getFileSize, getFileType, getFilename, getFilenameExtension, getLastAccessedOn, getLastModifiedOn, getOutputExtension, getPath, getPathname, gunzip, gzip, hasBeenModified, is.na, isFile, isGzipped, linkTo, readChecksum, renameTo, renameToUpperCaseExt, setAttribute, setAttributes, setAttributesBy, setAttributesByTags, setExtensionPattern, testAttributes, validate, validateChecksum, writeChecksum, directoryItem, directoryItems, directoryStructure, directoryStructure<-, getDataSetName, getOrganismName, getSampleName, loadCacheFile, loadCacheFileItem, memoizedCall2, saveCacheFile, saveCacheFileItem, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\section{Compression}{
 Currently, the package only supports non-compressed FASTA files.
}

\section{Filenames}{
 Currently, FASTA files with commas in their filenames should be avoided
 because they are not supported by Bowtie2.
}

\section{On different FASTA files}{
 Some FASTA files represent nucleotides as
 upper-case letters A, C, G, T (and N), whereas others as
 lower-case letters a, c, g, t (and n).
 The format of the sequence names differ between FASTA file
 (even for the same genome/organism), e.g. 'chr1' versus '1'.
}

\author{Henrik Bengtsson}

\seealso{
  [1] \url{http://www.wikipedia.org/wiki/FASTA_format}
}
\keyword{classes}
