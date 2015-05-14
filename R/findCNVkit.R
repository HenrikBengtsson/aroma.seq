findCNVkit <- function(..., command="cnvkit.py") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'command':
  command <- match.arg(command)

  # Aroma-specific variable
  path <- getExternalHome("CNVKIT_HOME")

  versionPattern <- c("version"="([0-9.]+)")
  res <- findExternal(command=command, path=path, versionPattern=versionPattern, ...)

  res
} # findCNVkit()

############################################################################
# HISTORY:
# 2015-05-13
# o Created.
############################################################################
