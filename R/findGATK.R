findGATK <- function(mustExist=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'mustExist':
  mustExist <- Arguments$getLogical(mustExist);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating GATK software");

  command <- "GATK";
  verbose && cat(verbose, "Command: ", command);

  # Check for cached results
  res <- .findCache(name=command);
  if (!is.null(res)) {
    pathname <- res$path;
    verbose && cat(verbose, "Found cached result.");
    verbose && exit(verbose);
    return(pathname);
  }

  path <- Sys.getenv("GATK_HOME");
  verbose && printf(verbose, "System variable 'GATK_HOME': '%s'\n", path);
  if (path == "") path <- NULL;

  if (is.null(path) || !isDirectory(path)) {
    if (mustExist) {
      throw(sprintf("Failed to located GATK home directory"));
    }
    return(NULL);
  }

  path <- Arguments$getReadablePath(path, mustWork=FALSE);
  verbose && cat(verbose, "Located directory: ", path);

  # Get main GATK jar file
  filename <- "GenomeAnalysisTK.jar";
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
  verbose && cat(verbose, "Located main jar file: ", pathname);


  if (isFile(pathname)) {
    .findCache(name=command, path=pathname);
    ver <- systemGATK("--version", stdout=TRUE, stderr=TRUE);
    attr(pathname, "version") <- ver;
    .findCache(name=command, path=pathname);
  } else if (mustExist) {
    throw(sprintf("Failed to located GATK tools"));
  }

  verbose && exit(verbose);

  pathname;
} # findGATK()
