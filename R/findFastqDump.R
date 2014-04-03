###########################################################################/**
# @RdocFunction findFastqDump
#
# @title "Locates the fastq-dump executable"
#
# \description{
#  @get "title" on the current system.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
#   \item{commandName}{'fastq-dump'; command name to find}
#   \item{versionPattern}{regexp to use if version not found properly; default should work}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
findFastqDump <- function(...,
                          commandName='fastq-dump',
                          versionPattern=".*fastq-dump[ ]*:[ ]*([0-9.]+)",
                          verbose=FALSE)
{
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  verbose && enter(verbose, "Locating FastqDump software");
  res <- findExternal(command=commandName, versionPattern=versionPattern, ...);
  res
} # findFastqDump()

############################################################################
# HISTORY:
# 2014-03-06
# o Copied from findFastQC(), findPerl()
############################################################################
