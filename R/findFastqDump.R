###########################################################################/**
# @RdocFunction findFastqDump
#
# @title "Locates the FastqDump executable"
#
# \description{
#  @get "title" on the current system.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
findFastqDump <- function(commandName='fastq-dump',...) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  verbose && enter(verbose, "Locating FastqDump software");
  versionPattern <- "fastq-dump[ ]*:[s ]*([0-9.]+)";
  findExternal(command=commandName, versionPattern=versionPattern, ...);
} # findFastqDump()

############################################################################
# HISTORY:
# 2014-03-06
# o Copied from findFastQC(), findPerl()
############################################################################
