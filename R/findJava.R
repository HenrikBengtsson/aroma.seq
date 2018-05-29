findJava <- function(...) {
  versionPattern <- c("-version"='.*version ["]?([0-9.]+)["]?.*')
  findExternal(command="java", versionPattern=versionPattern, ...)
} # findJava()
