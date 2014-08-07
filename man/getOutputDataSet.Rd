%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  getOutputDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getOutputDataSet}
\alias{getOutputDataSet}

\alias{getOutputDataSet.AromaSeqTransform}
\alias{getOutputDataSet.AbstractAlignment}
\alias{getOutputDataSet.BamDownsampler}
\alias{getOutputDataSet.Bowtie2Alignment}
\alias{getOutputDataSet.FastqDownsampler}
\alias{getOutputDataSet.PicardDuplicateRemoval}
\alias{getOutputDataSet.TopHat2Alignment}
\alias{getOutputDataSet.TotalCnBinnedCounting}

\title{Gets the (complete or incomplete) processed output data set}

\description{
  Gets the (complete or incomplete) processed output data set.
}

\usage{
 \method{getOutputDataSet}{AromaSeqTransform}(...)
 \method{getOutputDataSet}{AbstractAlignment}(this, onMissing=c("drop", "NA", "error"), ...)
 \method{getOutputDataSet}{Bowtie2Alignment}(this, onMissing=c("drop", "NA", "error"), ...)
 \method{getOutputDataSet}{FastqDownsampler}(this, ...)
 \method{getOutputDataSet}{PicardDuplicateRemoval}(this, onMissing=c("drop", "NA", "error"), ...)
 \method{getOutputDataSet}{TopHat2Alignment}(this, onMissing=c("drop", "NA", "error"), ...)
 \method{getOutputDataSet}{TotalCnBinnedCounting}(this, onMissing=c("drop", "NA", "error"), ...)
}

\arguments{
 \item{onMissing}{A \code{\link[base]{character}} string specifying how non-processed files
  should be returned.
  If \code{"drop"}, they are ignored and not part of the returned
  data set.
  If \code{"dropall"}, \code{\link[base]{NULL}} is returned unless all files are processed.
  If \code{"NA"}, they are represented as a "missing" file.
  If \code{"error"}, they are not accepted and an exception is thrown.
 }
 \item{...}{Not used.}
}

\value{
  Returns the output data set containing the same number of files as
  the input data set, except in the case where argument \code{onMissing}
  is \code{"drop"} or \code{"dropall"} and one or more files is not
  processed.
}

\seealso{
  This method is utilized by \code{\link{findFilesTodo}}().
}

\author{Henrik Bengtsson}


\keyword{methods}
\keyword{internal}
