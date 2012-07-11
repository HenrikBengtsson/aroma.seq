setMethodS3("bowtie2", "default", function(..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate external software
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Locating external software");
  command <- "bowtie2";
  verbose && cat(verbose, "Command: ", command)

  bin <- Sys.which(command);
  verbose && cat(verbose, "Located pathname: ", bin)

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call external software
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling external software");

  cmd <- sprintf("%s", shQuote(bin));
  verbose && cat(verbose, "System call: ", cmd)

  res <- system(cmd, ...);

  verbose && exit(verbose);

  invisible(res);
}) # bowtie2()


############################################################################
# HISTORY:
# 2012-07-11
# o Created bowtie2() stub.
############################################################################
