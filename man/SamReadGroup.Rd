%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SamReadGroup.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{SamReadGroup}
\docType{class}
\alias{SamReadGroup}


\title{The SamReadGroup class}

\description{
 Package:  aroma.seq \cr
\bold{Class SamReadGroup}\cr

\code{logical}\cr
\code{~~|}\cr
\code{~~+--}\code{\link[R.oo]{BasicObject}}\cr
\code{~~~~~~~|}\cr
\code{~~~~~~~+--}\code{SamReadGroup}\cr

\bold{Directly known subclasses:}\cr
\cr

public static class \bold{SamReadGroup}\cr
extends \link[R.oo]{BasicObject}\cr



 A SamReadGroup object represents a SAM Read Group.
}

\usage{
SamReadGroup(id=NULL, sample=NULL, library=NULL, platform=NULL, platformUnit=NULL,
  sequencingCenter=NULL, description=NULL, runDate=NULL, flowOrder=NULL,
  keySequence=NULL, program=NULL, predictedInsertSize=NULL, ...)
}

\arguments{
 \item{id}{Specifies the \code{ID} read group field.}
 \item{sample}{Specifies the \code{SM} read group field.}
 \item{library}{Specifies the \code{LB} read group field.}
 \item{platform}{Specifies the \code{PL} read group field.}
 \item{platformUnit}{Specifies the \code{PU} read group field.}
 \item{sequencingCenter}{Specifies the \code{CN} read group field.}
 \item{description}{Specifies the \code{DS} read group field.}
 \item{runDate}{Specifies the \code{DT} read group field.}
 \item{flowOrder}{Specifies the \code{FO} read group field.}
 \item{keySequence}{Specifies the \code{KS} read group field.}
 \item{program}{Specifies the \code{PG} read group field.}
 \item{predictedInsertSize}{Specifies the \code{PI} read group field.}
 \item{...}{Additional named arguments, including two-letter read
   group keys for the above, e.g. 'SM'.}
}

\section{Fields and Methods}{
 \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{asSamList} \tab  -\cr
 \tab \code{asString} \tab  -\cr
 \tab \code{byScanBamHeader} \tab  -\cr
 \tab \code{hasDescription} \tab  -\cr
 \tab \code{hasID} \tab  -\cr
 \tab \code{hasLibrary} \tab  -\cr
 \tab \code{hasPlatform} \tab  -\cr
 \tab \code{hasPlatformUnit} \tab  -\cr
 \tab \code{hasRunDate} \tab  -\cr
 \tab \code{hasSample} \tab  -\cr
 \tab \code{hasSequencingCenter} \tab  -\cr
 \tab \code{isEmpty} \tab  -\cr
 \tab \code{merge} \tab  -\cr
 \tab \code{validate} \tab  -\cr
}


 \bold{Methods inherited from BasicObject}:\cr
$, $<-, [[, [[<-, as.character, attach, detach, equals, extend, getFields, getInstantiationTime, hasField, hashCode, isReferable, objectSize, print, asThis

 \bold{Methods inherited from logical}:\cr
%*%,CsparseMatrix,numLike-method, %*%,Matrix,numLike-method, %*%,numLike,CsparseMatrix-method, %*%,numLike,Matrix-method, %*%,numLike,sparseVector-method, %*%,sparseVector,numLike-method, Arith,abIndex,numLike-method, Arith,ddenseMatrix,logical-method, Arith,ddiMatrix,logical-method, Arith,dgCMatrix,logical-method, Arith,dgeMatrix,logical-method, Arith,dpoMatrix,logical-method, Arith,dppMatrix,logical-method, Arith,dsparseMatrix,logical-method, Arith,lMatrix,logical-method, Arith,ldiMatrix,logical-method, Arith,logical,ddenseMatrix-method, Arith,logical,ddiMatrix-method, Arith,logical,dgCMatrix-method, Arith,logical,dgeMatrix-method, Arith,logical,dpoMatrix-method, Arith,logical,dppMatrix-method, Arith,logical,dsparseMatrix-method, Arith,logical,lMatrix-method, Arith,logical,ldiMatrix-method, Arith,logical,nMatrix-method, Arith,nMatrix,logical-method, Arith,numLike,abIndex-method, Compare,dMatrix,logical-method, Compare,lMatrix,logical-method, Compare,logical,dMatrix-method, Compare,logical,lMatrix-method, Compare,logical,nMatrix-method, Compare,nMatrix,logical-method, Logic,dMatrix,logical-method, Logic,lMatrix,logical-method, Logic,logical,dMatrix-method, Logic,logical,lMatrix-method, Logic,logical,nMatrix-method, Logic,nMatrix,logical-method, NSBS,logical-method, Ops,AtomicList,atomic-method, Ops,CompressedAtomicList,atomic-method, Ops,DelayedArray,vector-method, Ops,Rle,vector-method, Ops,SimpleAtomicList,atomic-method, Ops,atomic,AtomicList-method, Ops,atomic,CompressedAtomicList-method, Ops,atomic,SimpleAtomicList-method, Ops,atomicVector,sparseVector-method, Ops,ddiMatrix,logical-method, Ops,dpoMatrix,logical-method, Ops,dppMatrix,logical-method, Ops,ldiMatrix,logical-method, Ops,logical,dpoMatrix-method, Ops,logical,dppMatrix-method, Ops,nonStructure,vector-method, Ops,sparseVector,atomicVector-method, Ops,structure,vector-method, Ops,vector,DelayedArray-method, Ops,vector,Rle-method, Ops,vector,nonStructure-method, Ops,vector,structure-method, [,CsparseMatrix,index,index,logical-method, [,CsparseMatrix,index,missing,logical-method, [,CsparseMatrix,missing,index,logical-method, [,Matrix,index,index,missing-method, [,Matrix,logical,missing,missing-method, [,Matrix,missing,index,missing-method, [,TsparseMatrix,index,index,logical-method, [,TsparseMatrix,index,missing,logical-method, [,TsparseMatrix,missing,index,logical-method, [,abIndex,index,ANY,ANY-method, [,denseMatrix,index,index,logical-method, [,denseMatrix,index,missing,logical-method, [,denseMatrix,missing,index,logical-method, [,diagonalMatrix,index,index,logical-method, [,diagonalMatrix,index,missing,logical-method, [,diagonalMatrix,missing,index,logical-method, [,indMatrix,index,missing,logical-method, [,sparseMatrix,index,index,logical-method, [,sparseMatrix,index,missing,logical-method, [,sparseMatrix,missing,index,logical-method, [,sparseVector,index,ANY,ANY-method, [<-,CsparseMatrix,Matrix,missing,replValue-method, [<-,CsparseMatrix,index,index,replValue-method, [<-,CsparseMatrix,index,index,sparseVector-method, [<-,CsparseMatrix,index,missing,replValue-method, [<-,CsparseMatrix,index,missing,sparseVector-method, [<-,CsparseMatrix,matrix,missing,replValue-method, [<-,CsparseMatrix,missing,index,replValue-method, [<-,CsparseMatrix,missing,index,sparseVector-method, [<-,Matrix,ldenseMatrix,missing,replValue-method, [<-,Matrix,lsparseMatrix,missing,replValue-method, [<-,Matrix,matrix,missing,replValue-method, [<-,Matrix,ndenseMatrix,missing,replValue-method, [<-,Matrix,nsparseMatrix,missing,replValue-method, [<-,RsparseMatrix,index,index,replValue-method, [<-,RsparseMatrix,index,index,sparseVector-method, [<-,RsparseMatrix,index,missing,replValue-method, [<-,RsparseMatrix,index,missing,sparseVector-method, [<-,RsparseMatrix,matrix,missing,replValue-method, [<-,RsparseMatrix,missing,index,replValue-method, [<-,RsparseMatrix,missing,index,sparseVector-method, [<-,TsparseMatrix,Matrix,missing,replValue-method, [<-,TsparseMatrix,index,index,replValue-method, [<-,TsparseMatrix,index,index,sparseVector-method, [<-,TsparseMatrix,index,missing,replValue-method, [<-,TsparseMatrix,index,missing,sparseVector-method, [<-,TsparseMatrix,matrix,missing,replValue-method, [<-,TsparseMatrix,missing,index,replValue-method, [<-,TsparseMatrix,missing,index,sparseVector-method, [<-,denseMatrix,index,index,replValue-method, [<-,denseMatrix,index,missing,replValue-method, [<-,denseMatrix,matrix,missing,replValue-method, [<-,denseMatrix,missing,index,replValue-method, [<-,diagonalMatrix,index,index,replValue-method, [<-,diagonalMatrix,index,index,sparseMatrix-method, [<-,diagonalMatrix,index,index,sparseVector-method, [<-,diagonalMatrix,index,missing,replValue-method, [<-,diagonalMatrix,index,missing,sparseMatrix-method, [<-,diagonalMatrix,index,missing,sparseVector-method, [<-,diagonalMatrix,matrix,missing,replValue-method, [<-,diagonalMatrix,missing,index,replValue-method, [<-,diagonalMatrix,missing,index,sparseMatrix-method, [<-,diagonalMatrix,missing,index,sparseVector-method, [<-,indMatrix,index,ANY,ANY-method, [<-,indMatrix,missing,index,ANY-method, [<-,sparseVector,index,missing,replValueSp-method, [<-,sparseVector,sparseVector,missing,replValueSp-method, append,Rle,vector-method, append,vector,Rle-method, as.data.frame, as.raster, bpexportglobals<-,BiocParallelParam,logical-method, bplog<-,SerialParam,logical-method, bplog<-,SnowParam,logical-method, bpprogressbar<-,BiocParallelParam,logical-method, bpstopOnError<-,BiocParallelParam,logical-method, cbind2,Matrix,atomicVector-method, cbind2,atomicVector,Matrix-method, cbind2,atomicVector,ddiMatrix-method, cbind2,atomicVector,ldiMatrix-method, cbind2,ddiMatrix,atomicVector-method, cbind2,ldiMatrix,atomicVector-method, coerce,ANY,logical-method, coerce,Matrix,logical-method, coerce,Rle,logical-method, coerce,Vector,logical-method, coerce,abIndex,vector-method, coerce,atomicVector,dsparseVector-method, coerce,atomicVector,sparseVector-method, coerce,dgCMatrix,vector-method, coerce,dsCMatrix,vector-method, coerce,logical,IRanges-method, coerce,logical,LLint-method, coerce,logical,NormalIRanges-method, coerce,logical,abIndex-method, coerce,numLike,dgeMatrix-method, coerce,sparseVector,logical-method, coerce,vector,Grouping-method, coerce,vector,ManyToManyGrouping-method, coerce,vector,ManyToOneGrouping-method, crossprod,CsparseMatrix,numLike-method, crossprod,Matrix,numLike-method, crossprod,dgeMatrix,numLike-method, crossprod,numLike,CsparseMatrix-method, crossprod,numLike,Matrix-method, crossprod,numLike,dgeMatrix-method, crossprod,numLike,sparseVector-method, crossprod,sparseVector,numLike-method, elementType,vector-method, extractROWS,vector_OR_factor,RangeNSBS-method, extractROWS,vector_OR_factor,RangesNSBS-method, invertStrand,logical-method, match,CompressedList,vector-method, match,List,vector-method, match,XStringSet,vector-method, match,vector,XStringSet-method, mstack,vector-method, pcompare,XStringSet,vector-method, pcompare,vector,XStringSet-method, pmax2,DelayedArray,vector-method, pmax2,vector,DelayedArray-method, pmin2,DelayedArray,vector-method, pmin2,vector,DelayedArray-method, rbind2,Matrix,atomicVector-method, rbind2,atomicVector,Matrix-method, rbind2,atomicVector,ddiMatrix-method, rbind2,atomicVector,ldiMatrix-method, rbind2,ddiMatrix,atomicVector-method, rbind2,ldiMatrix,atomicVector-method, rename,vector-method, shiftApply,vector,vector-method, splitRanges,vector_OR_factor-method, strand,logical-method, tcrossprod,CsparseMatrix,numLike-method, tcrossprod,Matrix,numLike-method, tcrossprod,dgeMatrix,numLike-method, tcrossprod,numLike,CsparseMatrix-method, tcrossprod,numLike,Matrix-method, tcrossprod,numLike,dgeMatrix-method, tcrossprod,numLike,sparseVector-method, tcrossprod,sparseVector,numLike-method, whichVector, window<-,vector-method


}

\author{Henrik Bengtsson}

\references{
 [1] The SAM Format Specification Working Group,
     \emph{The SAM Format Specification}, Sept 7, 2011.\cr
 [2] \emph{Collected FAQs about BAM files: What's the meaning
           of the standard read group fields?}, GATK Forums,
     Mar 2013. \url{http://goo.gl/PWv8MX} \cr
}
\keyword{classes}
