%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Do not modify this file since it was automatically generated from:
% 
%  FastaReferenceFile.R
% 
% by the Rdoc compiler part of the R.oo package.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\name{buildBwaIndexSet.FastaReferenceFile}
\alias{buildBwaIndexSet.FastaReferenceFile}
\alias{FastaReferenceFile.buildBwaIndexSet}
\alias{buildBwaIndexSet,FastaReferenceFile-method}

\title{Builds a BWA index files set}

\description{
  Builds a BWA index files set.
}

\usage{
\method{buildBwaIndexSet}{FastaReferenceFile}(this, method, ..., skip=TRUE, verbose=FALSE)
}

\arguments{
 \item{method}{A \code{\link[base]{character}} string specifying the algorithm to use,
    cf. \code{\link{bwaIndexPrefix}}()}
 \item{...}{Additional arguments passed to \code{\link{bwaIndex}}().}
 \item{skip}{If \code{\link[base:logical]{TRUE}}, the index files are not rebuilt if already available.}
 \item{verbose}{See \code{\link[R.utils]{Verbose}}.}
}

\value{
  Returns a \code{\link[R.filesets]{GenericDataFileSet}} consisting of the BWA index files.
}

\seealso{
  Internally, \code{\link{bwaIndex}}() is used.
}

\references{
 [1] Edwards Bioinformatics Lab,
     \emph{How to create a database for BWA and BWA-SW}, 2013.
     \url{http://edwards.sdsu.edu/research/index.php/robert/282-how-to-create-a-database-for-bwa-and-bwa-sw} \cr
 [2] Henrik Bengtsson, \emph{bwa index -a is: Details on database 2GB limit?},
     bwa-help thread on 2013-11-18.
     \url{https://sourceforge.net/mailarchive/message.php?msg_id=31649355}\cr
}

\author{Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
