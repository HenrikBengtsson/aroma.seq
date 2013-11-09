%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastqDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{FastqDataSet}
\docType{class}
\alias{FastqDataSet}


\title{The FastqDataSet class}

\description{
 Package:  aroma.seq \cr
\bold{Class FastqDataSet}\cr

\code{\link[R.oo]{Object}}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.filesets]{FullNameInterface}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{\link[R.filesets]{GenericDataFileSet}}\cr
\code{~~~~~~~~~~~~|}\cr
\code{~~~~~~~~~~~~+--}\code{FastqDataSet}\cr

\bold{Directly known subclasses:}\cr
\link[aroma.seq]{IlluminaFastqDataSet}, \link[aroma.seq]{PairedEndFastqDataSet}\cr

public static class \bold{FastqDataSet}\cr
extends \link[R.filesets]{GenericDataFileSet}\cr



 An FastqDataSet object represents a set of \code{\link{FastqDataFile}}:s.
}

\usage{
FastqDataSet(files=NULL, ..., paired=FALSE)
}

\arguments{
  \item{files}{A \code{\link[base]{list}} of \code{\link{FastqDataFile}}:s.}
  \item{...}{Arguments passed to \code{\link[R.filesets]{GenericDataFileSet}}.}
  \item{paired}{If \code{\link[base:logical]{TRUE}}, the data set contains paired-end reads,
    otherwise not.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{byName} \tab  -\cr
 \tab \code{doBWA} \tab  -\cr
 \tab \code{doBowtie2} \tab  -\cr
 \tab \code{doQDNAseq} \tab  -\cr
 \tab \code{doTopHat2} \tab  -\cr
 \tab \code{findByName} \tab  -\cr
 \tab \code{getDefaultSamReadGroup} \tab  -\cr
 \tab \code{getSamReadGroup} \tab  -\cr
 \tab \code{report} \tab  -\cr
 \tab \code{setSamReadGroup} \tab  -\cr
}


 \bold{Methods inherited from GenericDataFileSet}:\cr
[, [[, anyNA, append, appendFiles, appendFullNamesTranslator, appendFullNamesTranslatorByNULL, appendFullNamesTranslatorByTabularTextFile, appendFullNamesTranslatorByTabularTextFileSet, appendFullNamesTranslatorBydata.frame, appendFullNamesTranslatorByfunction, appendFullNamesTranslatorBylist, as.character, as.list, byName, byPath, c, clearCache, clearFullNamesTranslator, clone, copyTo, equals, extract, findByName, getChecksum, getChecksumObjects, getDefaultFullName, getFile, getFileClass, getFileSize, getFiles, getFullNames, getNames, getOneFile, getPath, getPathnames, getSubdirs, gunzip, gzip, hasFile, indexOf, nbrOfFiles, resetFullNames, seq, setFullNamesTranslator, sortBy, update2, updateFullName, updateFullNames, validate, dsApply, getFullNameTranslatorSet, getParentName

 \bold{Methods inherited from FullNameInterface}:\cr
appendFullNameTranslator, appendFullNameTranslatorByNULL, appendFullNameTranslatorByTabularTextFile, appendFullNameTranslatorByTabularTextFileSet, appendFullNameTranslatorBycharacter, appendFullNameTranslatorBydata.frame, appendFullNameTranslatorByfunction, appendFullNameTranslatorBylist, clearFullNameTranslator, clearListOfFullNameTranslators, getDefaultFullName, getFullName, getFullNameTranslator, getListOfFullNameTranslators, getName, getTags, hasTag, hasTags, resetFullName, setFullName, setFullNameTranslator, setListOfFullNameTranslators, setName, setTags, updateFullName

 \bold{Methods inherited from Object}:\cr
$, $<-, [[, [[<-, as.character, attach, attachLocally, clearCache, clearLookupCache, clone, detach, equals, extend, finalize, gc, getEnvironment, getFieldModifier, getFieldModifiers, getFields, getInstantiationTime, getStaticInstance, hasField, hashCode, ll, load, objectSize, print, save, findRspReportTemplate, asThis


}

\section{Paired reads}{
  There is a community/industry convention for paired-end runs:
  "The reads are reported two FASTQ files, such that the n:th read in
  the first file is mate-paired to the n:th read in the second file.
  The read IDs must match." [1]

  If argument \code{paired=TRUE}, the \code{files} arguments is assumed
  to contain the "R1" files whereas the corresponding "R2" files are
  implicit (inferred and located by matching the filenames).
}

\references{
 [1] Simon Anders,
     \emph{High-throughput sequencing: Alignment and related topic},
     (38 slides), EMBL Heidelberg, June 2013.\cr
}

\author{Henrik Bengtsson}
\keyword{classes}
