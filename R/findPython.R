findPython <- function(...) {
  versionPattern <- c("--version"='.* ([0-9.a-z]+)?');
  findExternal(command="python", versionPattern=versionPattern, ...);
} # findPython()


############################################################################
# HISTORY:
# 2013-10-30
# o Created from findPerl.R.
############################################################################
