%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  getOutputDataSet.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{getOutputDataSet}
\alias{getOutputDataSet}

\alias{getOutputDataSet.AbstractAlignment}
\alias{getOutputDataSet.FastqDownsampler}
\alias{getOutputDataSet.PicardDuplicateRemoval}
\alias{getOutputDataSet.TopHat2Alignment}
\alias{getOutputDataSet.TotalCnBinnedCounting}

\title{Gets the (complete or incomplete) processed output data set}

\description{
  Gets the (complete or incomplete) processed output data set.
}

\usage{
getOutputDataSet(...)
}

\arguments{
 \item{onMissing}{A \code{\link[base]{character}} string specifying how non-processed files
  should be returned.
  If \code{"drop"}, they are ignored and not part of the returned
  data set.
  If \code{"NA"}, they are represented as a "missing" file.
  If \code{"error"}, they are not accepted and an exception is thrown.
 }
 \item{...}{Not used.}
}

\value{
  Returns the output data set containing the same number of files as
  the input data set, except in the case where \code{onMissing="drop"}.
}

\seealso{
  This method is utilized by \code{\link{findFilesTodo}}().
}

\author{Henrik Bengtsson}


\keyword{methods}
\keyword{internal}
