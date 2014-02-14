%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  Bowtie2IndexSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{Bowtie2IndexSet}
\docType{class}
\alias{Bowtie2IndexSet}


\title{The Bowtie2IndexSet class}

\description{
 Package:  aroma.seq \cr
\bold{Class Bowtie2IndexSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{\link[aroma.seq]{AbstractIndexSet}}\cr
\code{~~~~~~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~~~~~~+--}\code{Bowtie2IndexSet}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{Bowtie2IndexSet}\cr
extends \emph{\link[aroma.seq]{AbstractIndexSet}}\cr



 An Bowtie2IndexSet object represents a set of \code{\link{Bowtie2IndexFile}}:s.
}

\usage{
Bowtie2IndexSet(...)
}

\arguments{
  \item{...}{Arguments passed to \code{\link{AbstractIndexSet}}.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{\link[aroma.seq:buildTopHat2TranscriptomeIndexSet.Bowtie2IndexSet]{buildTopHat2TranscriptomeIndexSet}} \tab Calls the TopHat executable to build a transcriptome index on a gene model in a reference genome.\cr
 \tab \code{getSequenceNames} \tab  -\cr
 \tab \code{getSummary} \tab  -\cr
 \tab \code{isComplete} \tab  -\cr
}


 \bold{Methods inherited from AbstractIndexSet}:\cr
as.character, byPrefix, getFastaReferenceFile, getIndexPrefix, getOrganism, isComplete

 \bold{Methods inherited from GenericDataFileSet}:\cr
[, [[, anyDuplicated, anyNA, append, appendFiles, appendFullNamesTranslator, appendFullNamesTranslatorByNULL, appendFullNamesTranslatorByTabularTextFile, appendFullNamesTranslatorByTabularTextFileSet, appendFullNamesTranslatorBydata.frame, appendFullNamesTranslatorByfunction, appendFullNamesTranslatorBylist, as.character, as.list, byName, byPath, c, clearCache, clearFullNamesTranslator, clone, copyTo, duplicated, equals, extract, findByName, findDuplicated, getChecksum, getChecksumFileSet, getChecksumObjects, getDefaultFullName, getFile, getFileClass, getFileSize, getFiles, getFullNames, getNames, getOneFile, getPath, getPathnames, getSubdirs, gunzip, gzip, hasFile, indexOf, is.na, na.omit, nbrOfFiles, resetFullNames, seq, setFullNamesTranslator, sortBy, unique, update2, updateFullName, updateFullNames, validate, directoryItem, directoryItems, directoryStructure, directoryStructure<-, dsApply, updateDirectoryStructure, getFullNameTranslatorSet, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\author{Henrik Bengtsson}

\references{
 [1] Trapnell et al. \emph{Differential gene and transcript expression
     analysis of RNA-seq experiments with TopHat and Cufflinks}.
     Nat Protoc, 2012.\cr
}


\keyword{classes}
\keyword{internal}
