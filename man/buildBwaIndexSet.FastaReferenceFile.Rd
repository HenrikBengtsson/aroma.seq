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
\method{buildBwaIndexSet}{FastaReferenceFile}(this, method=c("bwtsw", "is"), ..., skip=TRUE, verbose=FALSE)
}

\arguments{
 \item{method}{A \code{\link[base]{character}} string specifying the algorithm to use for
    building the index set. All methods gives identical results [1].
    The default is such that it can handle also large genomes, including
    the human genome.}
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
 [1] Thread \emph{bwa index option bwtsw}, SEQanswers, 2010-07-13.
     \url{http://seqanswers.com/forums/showthread.php?t=5921}\cr
 [2] Edwards Bioinformatics Lab,
     \emph{How to create a database for BWA and BWA-SW}, 2013.
     \url{http://edwards.sdsu.edu/research/index.php/robert/282-how-to-create-a-database-for-bwa-and-bwa-sw} \cr
 [3] Henrik Bengtsson, \emph{bwa index -a is: Details on database 2GB limit?},
     bwa-help thread on 2013-11-18.
     \url{https://sourceforge.net/mailarchive/message.php?msg_id=31649355}\cr
}

\author{Henrik Bengtsson}
\keyword{internal}
\keyword{methods}
