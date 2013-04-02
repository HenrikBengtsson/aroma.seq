findBWA <- function(...) {
  versionPattern <- c("Version:[ ]*([0-9.-_]+).*");
  findExternal(command="bwa", versionPattern=versionPattern, ...);
}


############################################################################
# HISTORY:
# 2013-04-01
# o Now findBWA() sets attribute 'version', iff possible.
# 2012-09-24
# o Created.
############################################################################
