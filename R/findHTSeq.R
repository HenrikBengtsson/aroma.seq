###########################################################################/**
# @RdocFunction findHTSeq
#
# @title "Locates the HTSeq executable"
#
# \description{
#  @get "title" on the current system.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "findExternal".}
#   \item{command}{The type of HTSeq executable to locate.}
# }
#
# \details{
#  The HTSeq tool is search for as follows:
#  \enumerate{
#   \item ...
#  }
# }
#
# @author "HB,TT"
#
# \seealso{
#  [1] HTSeq: Analysing high-throughput sequencing data with Python,
#      June 2013.
#      \url{http://www-huber.embl.de/users/anders/HTSeq/}
# }
#
# @keyword internal
#*/###########################################################################
findHTSeq <- function(..., command=c("htseq-count", "htseq-qa")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'command':
  command <- match.arg(command);

  versionPattern <- c("-version"=".*version ([0-9.]+).*");
  findExternal(command=command, versionPattern=versionPattern, ...);
} # findHTSeq()

############################################################################
# HISTORY:
# 2013-06-20 [HB]
# o Added Rdoc help.
# o Renamed from findHtseq() to findHTSeq().
# 2013-05-31
# o TAT:  Created findHtseq
############################################################################


## Punt on version checking for now.  Need to brush up on Python first.
## Version info is avail e.g. via
## ./HTSeq/_version.py:1:__version__ = "0.5.4p3"
## also via the VERSION file in the install dir.
## Something like this might work:
##  pkgVersion <- sub("[^0-9.].*", "", system("python -c 'import HTSeq; print HTSeq.__version__'", intern=TRUE))

