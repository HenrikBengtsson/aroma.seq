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

