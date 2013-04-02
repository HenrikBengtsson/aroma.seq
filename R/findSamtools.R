findSamtools <- function(...) {
  versionPattern <- c("Version:[ ]*([0-9.-]+).*");
  findExternal(command="samtools", versionPattern=versionPattern, ...);
} # findSamtools()


############################################################################
# HISTORY:
# 2012-09-25
# o Created.
############################################################################
