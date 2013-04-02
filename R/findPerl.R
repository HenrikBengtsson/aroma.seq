findPerl <- function(...) {
  versionPattern <- c("-version"='.*This is .* [(]?v([0-9.]+)[)]?.*');
  findExternal(command="perl", versionPattern=versionPattern, ...);
} # findPerl()


############################################################################
# HISTORY:
# 2013-04-01
# o Now findPerl() sets attribute 'version', iff possible.
# 2012-09-28
# o Created.
############################################################################
