findPython <- function(...) {
  versionPattern <- c("--version"='.* ([0-9.a-z]+)?');
  findExternal(command="python", versionPattern=versionPattern, ...);
} # findPython()
