findHTSeq <- function(..., command=c("htseq-count", "htseq-qa")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'command':
  command <- match.arg(command);

  versionPattern <- c(".*version ([0-9.]+(|p[0-9]+)).*");
  findExternal(command=command, versionPattern=versionPattern, ...);
} # findHTSeq()

############################################################################
# HISTORY:
# 2014-01-24 [HB]
# o BUG FIX: findHTSeq() failed to identify the version.
# 2013-06-20 [HB]
# o Renamed from findHtseq() to findHTSeq().
# 2013-05-31 [TT]
# o Created.
############################################################################
