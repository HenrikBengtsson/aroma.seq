findBWA <- function(...) {
  # Aroma-specific variable
  path <- getExternalHome("BWA_HOME");
  versionPattern <- c("Version:[ ]*([0-9.-_]+).*");
  findExternal(command="bwa", path=path, versionPattern=versionPattern, ...);
}
