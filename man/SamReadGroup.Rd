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

 \usage{SamReadGroup(sample=NULL, library=NULL, platform=NULL, platformUnit=NULL, sequencingCenter=NULL, description=NULL, runDate=NULL, flowOrder=NULL, keySequence=NULL, program=NULL, predictedInsertSize=NULL, ...)}

 \arguments{
  \item{sample}{Specifies the \code{SM} read group.}
  \item{library}{Specifies the \code{LB} read group.}
  \item{platform}{Specifies the \code{PL} read group.}
  \item{platformUnit}{Specifies the \code{PU} read group.}
  \item{sequencingCenter}{Specifies the \code{CN} read group.}
  \item{description}{Specifies the \code{DS} read group.}
  \item{runDate}{Specifies the \code{DT} read group.}
  \item{flowOrder}{Specifies the \code{FO} read group.}
  \item{keySequence}{Specifies the \code{KS} read group.}
  \item{program}{Specifies the \code{PG} read group.}
  \item{predictedInsertSize}{Specifies the \code{PI} read group.}
  \item{...}{Additional named arguments, including two-letter read
    group keys for the above, e.g. 'SM'.}
 }

 \section{Fields and Methods}{
  \bold{Methods:}\cr
\tabular{rll}{
 \tab \code{as.character} \tab  -\cr
 \tab \code{asBwaString} \tab  -\cr
 \tab \code{asSamList} \tab  -\cr
 \tab \code{asString} \tab  -\cr
 \tab \code{byScanBamHeader} \tab  -\cr
 \tab \code{hasDescription} \tab  -\cr
 \tab \code{hasLibrary} \tab  -\cr
 \tab \code{hasPlatform} \tab  -\cr
 \tab \code{hasPlatformUnit} \tab  -\cr
 \tab \code{hasRunDate} \tab  -\cr
 \tab \code{hasSample} \tab  -\cr
 \tab \code{hasSequencingCenter} \tab  -\cr
 \tab \code{merge} \tab  -\cr
 \tab \code{validate} \tab  -\cr
}


 \bold{Methods inherited from BasicObject}:\cr
asThis, $, $<-, [[, [[<-, as.character, attach, detach, equals, extend, getFields, getInstantiationTime, hasField, hashCode, isReferable, objectSize, print

 \bold{Methods inherited from logical}:\cr
as.data.frame, as.raster, whichVector


 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}

 \references{
  [1] The SAM Format Specification Working Group,
      \emph{The SAM Format Specification}, Sept 7, 2011.\cr
 }
\keyword{classes}
