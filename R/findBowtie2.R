findBowtie2 <- function(..., command=c("bowtie2", "bowtie2-align", "bowtie2-build", "bowtie2-inspect")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'command':
  command <- match.arg(command);

  versionPattern <- c("-version"=".*version ([0-9.]+).*");
  findExternal(command=command, versionPattern=versionPattern, ...);
} # findBowtie2()

############################################################################
# HISTORY:
# 2013-04-01
# o BUG FIX: findBowtie2() was incorrectly hardwired to 'bowtie2-align'.
# o Now findBowtie2() sets attribute 'version', iff possible.
# 2012-09-27
# o Now looking for bowtie2-align instead of bowtie2, because the latter
#   is a perl script calling the former.
# 2012-09-24
# o Created.
############################################################################
