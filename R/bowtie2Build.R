#
# Command line call:
#  bowtie2-build path/to/lambda_virus.fa lambda_virus
#
# R system call:
#  bowtie2Build("path/to/lambda_virus.fa");
#  bowtie2Build("path/to/lambda_virus.fa", outPath="path/to/", name="lambda_virus");
#  => calls =>
#  systemBowtie2Build("path/to/lambda_virus.fa", "lambda_virus");
#  => calls =>
#  system("bowtie2-build path/to/lambda_virus.fa lambda_virus");
#
setMethodS3("bowtie2Build", "default", function(inPathname, outPath=NULL, name=NULL, ..., overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'inPathname':
  inPathname <- Arguments$getReadablePathname(inPathname);

  # Argument 'outPath':
  if (is.null(outPath)) {
    outPath <- dirname(inPathname);   # e.g. /path/to
  }
  outPath <- Arguments$getWritablePath(outPath);

  # Argument 'name':
  if (is.null(name)) {
    inFilename <- basename(inPathname);   # e.g. lambda_virus.fa
    name <- gsub("[.][^.]$", ".bam", inFilename);  # e.g. lambda_virus
  }
  name <- Arguments$getCharacter(name);

  # Additional arguments
  args <- list(...);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running bowtie2-build");
  verbose && cat(verbose, "Input pathname: ", inPathname);
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, args);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate external software
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Locating external software");
  command <- "bowtie2-build";
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

  # Argument to bowtie2-index (maybe)
  outPrefix <- file.path(outPath, name); # e.g. path/to/lambda_virus

  cmdSwitches <- args$cmdSwitches
  binWithArgs <- paste(bin, cmdSwitches, inPathname, name)

  ## cmd <- sprintf("%s", shQuote(binWithArgs));
  cmd <- sprintf("%s", binWithArgs)  ## In emacs/ess, shQuote() version with single quotes added does not run

  verbose && cat(verbose, "System call: ", cmd)

  res <- system(cmd, ...);

  verbose && exit(verbose);

  verbose && exit(verbose);

  invisible(res);
}) # bowtie2-build()


############################################################################
# HISTORY:
# 2012-07-18
# o Created bowtie2-build() stub.
############################################################################
