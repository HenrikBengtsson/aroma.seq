findPerl <- function(...) {
  versionPattern <- c("-version"='.*This is .* [(]?v([0-9.]+)[)]?.*');
  findExternal(command="perl", versionPattern=versionPattern, ...);
} # findPerl()
