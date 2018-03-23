findSraToolkit <- function(command="fastq-dump", ...) {
  # Aroma-specific variable
  path <- getExternalHome("SRATOOLKIT_HOME");
  if (!is.null(path)) path <- file.path(path, "bin")

  if (is.null(path) || !isDirectory(path)) {
    pathnameT <- Sys.which(command)
    if (isFile(pathnameT)) {
      path <- dirname(pathnameT)
    }
  }

  versionPattern <- sprintf("%s[ ]*:[ ]*([0-9.-]+).*$", command);
  res <- findExternal(command=command, path=path, versionPattern=versionPattern, ...);
  res
}
