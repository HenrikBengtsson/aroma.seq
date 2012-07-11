setMethodS3("bowtie2", "default", function(inPathname, outPathname=NULL, ..., overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'inPathname':
  inPathname <- Arguments$getReadablePathname(inPathname);

  # Argument 'outPathname':
  if (is.null(outPathname)) {
    outPathname <- gsub("[.][^.]$", ".bam", inPathname);
  }
  outPathname <- Arguments$getWritablePathname(outPathname, mustNotExist=!overwrite);

  # Additional arguments
  args <- list(...);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running bowtie2");
  verbose && cat(verbose, "Input pathname: ", inPathname);
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, args);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate external software
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Locating external software");
  command <- "bowtie2";
  verbose && cat(verbose, "Command: ", command)

  bin <- Sys.which(command);
  verbose && cat(verbose, "Located pathname: ", bin)
  
  # Assert existence
  if (identical(bin, "") || !isFile(bin)) {
    throw("Failed to located external software (via the system PATH): ", command);
  }
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call external software
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling external software");

  cmd <- sprintf("%s", shQuote(bin));
  verbose && cat(verbose, "System call: ", cmd)

  res <- system(cmd, ...);

  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(res);
}) # bowtie2()


############################################################################
# HISTORY:
# 2012-07-11
# o Created bowtie2() stub.
############################################################################
