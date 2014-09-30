findSraToolkit <- function(command="fastq-dump", ...) {
  # Aroma-specific variable
  path <- getExternalHome("SRATOOLKIT_HOME");
  path <- file.path(path, "bin")
  versionPattern <- c("[ ]*([0-9.-]+).*");
  res <- findExternal(command=command, path=path, versionPattern=versionPattern, ...);
  res
}


############################################################################
# HISTORY:
# 2014-09-29
# o Created.
############################################################################
