findJava <- function(...) {
  versionPattern <- c("-version"='.*version ["]?([0-9.]+)["]?.*');
  findExternal(command="java", versionPattern=versionPattern, ...);
} # findJava()


############################################################################
# HISTORY:
# 2013-04-01
# o Now findJava() sets attribute 'version', iff possible.
# 2012-09-28
# o Created.
############################################################################
