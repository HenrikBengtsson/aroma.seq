findFastqDump <- function(...,
                          commandName='fastq-dump',
                          versionPattern=".*fastq-dump[ ]*:[ ]*([0-9.]+)",
                          verbose=FALSE)
{
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Aroma-specific variable
  path <- getExternalHome("FASTQDUMP_HOME");

  verbose && enter(verbose, "Locating FastqDump software");
  res <- findExternal(command=commandName, path=path, versionPattern=versionPattern, ...);
  res
} # findFastqDump()
