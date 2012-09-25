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

setMethodS3("bowtie2Build", "default", function(inPathnames,   ## vector of filenames
                                                outPath=NULL,  ## path to index
                                                outName=NULL,  ## basename of index
                                                optionsList,   ## list of bowtie2 options, for power users
                                                ...,
                                                overwrite=FALSE,
                                                verbose=FALSE,
                                                command="bowtie2-build") {

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Validate arguments
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  inPathnames <- sapply(inPathnames, Arguments$getReadablePathname)
  inPathsStr <- paste(inPathnames, collapse=",")   ## bowtie2 wants a comma-delimited string

  ## Default for outPath is the inPath
  if (is.null(outPath)) {
    outPath <- dirname(inPathnames[1])
  }
  outPath <- Arguments$getWritablePath(outPath);
  if (is.null(outPath)) {
    outName <- gsub("\\.fa$", "", basename(inPathnames[1]))   ## e.g. lambda_virus
    ## - Using the first infile may be confusing, if e.g. it is a subscripted name like "_1"; oh well
  }
  outPrefix <- file.path(outPath, outName); ## e.g. path/to/lambda_virus
  outPrefix <- Arguments$getWriteablePathname(outPrefix, mustNotExist=T)

  ## Additional arguments
  args <- list(...);  ## Not used here; just passed along

  ## Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  verbose && enter(verbose, paste("Running", command));
  verbose && cat(verbose, "Input pathnames: ", inPathsStr);
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, args);

  ## Locate external software
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Locating external software");
  verbose && cat(verbose, "Command: ", command)
  bin <- Sys.which(command);
  verbose && cat(verbose, "Located pathname: ", bin)

  ## Assert existence
  if (identical(bin, "") || !isFile(bin)) {
    throw("Failed to located external software (via the system PATH): ", command);
  }
  verbose && exit(verbose);

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ## Call external software
  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling external software");

  res <- systemBowtie2Build(inPathsStr, outPrefix, optionsList, ..., bin=bin)

  verbose && exit(verbose);
  verbose && exit(verbose);
  invisible(res);
}) ## bowtie2-build()


############################################################################
# HISTORY:
# 2012-07-18
# o Created bowtie2-build() stub. (HB)
# 2012-09-14
# o First real draft of "level 2" code (TAT)
############################################################################

