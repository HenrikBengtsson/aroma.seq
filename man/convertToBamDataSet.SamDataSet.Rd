%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  SamDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 \name{convertToBamDataSet.SamDataSet}
\alias{convertToBamDataSet.SamDataSet}
\alias{SamDataSet.convertToBamDataSet}
\alias{convertToBamDataSet.SamDataSet}
\alias{convertToBamDataSet,SamDataSet-method}


 \title{Converts the SAM set into a BAM set}

 \description{
   Converts the SAM set into a BAM set, where each BAM file is sorted and indexed.
 }

 \usage{\method{convertToBamDataSet}{SamDataSet}(this, path=getPath(this), ..., verbose=FALSE)}

 \arguments{
  \item{path}{The destination path.}
  \item{...}{Additional arguments passed to \code{convertToBamDataFile()}
   for each \code{\link{SamDataFile}}.}
  \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
 }

 \value{
   Returns a \code{\link{BamDataSet}}.
 }

 \seealso{
   Internally \code{\link[Rsamtools]{asBam}} is utilized.
 }

 \author{Henrik Bengtsson (\url{http://www.braju.com/R/})}
\keyword{internal}
\keyword{methods}
