###########################################################################/**
# @RdocGeneric gatk
# @alias gatk.default
#
# @title "Calls the GATK executable"
#
# \description{
#  @get "title".
# }
#
# \usage{
#  @usage gatk,default
# }
#
# \arguments{
#   \item{...}{(Not used)}.
#   \item{outPath}{Directory where result files are written.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "HB"
#
# \references{
#  [1] GATK: The Genome Analysis Toolkit,
#      Broad Institute, 2014.
#      \url{http://www.broadinstitute.org/gatk/}
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("gatk", "default", function(..., analysisType=NULL, inputFile=NULL, referenceSequence=NULL, interval=NULL, outPath="gatkData/", verbose=FALSE) {
  # Make sure to evaluate registered onExit() statements
  on.exit(eval(onExit()));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'inputFile'
  if (length(inputFile) > 0L) {
    inputFile <- Arguments$getReadablePathnames(inputFile, absolute=TRUE);
    assertNoDuplicated(inputFile);
  }

  # Argument 'outPath'
  outPath <- Arguments$getWritablePath(outPath, mustNotExist=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
     pushState(verbose);
     on.exit(popState(verbose), add=TRUE);
  }

  verbose && enter(verbose, "Running gatk()");
  verbose && cat(verbose, "Input files:");
  verbose && print(verbose, sapply(inputFile, FUN=getRelativePath));
  verbose && cat(verbose, "Output directory: ", outPath);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate output atomically by writing to a temporary directory
  # that is renamed upon completion.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  outPathOrg <- outPath;
  outPath <- sprintf("%s.TMP", outPath);
  outPath <- Arguments$getWritablePath(outPath, mustNotExist=TRUE);
  verbose && cat(verbose, "Temporary output directory: ", outPath);
  # At the end, assume failure, unless successful.
  outPathFinal <- sprintf("%s.ERROR", outPath);
  onExit({
    removeDirectory(outPathOrg, recursive=FALSE, mustExist=FALSE);
    file.rename(outPath, outPathFinal);
  });


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # While running, let the output directory be the working directory.
  # Make sure the working directory is restored when exiting.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opwd <- setwd(outPath);
  verbose && cat(verbose, "Original working directory: ", opwd);
  onExit({
    verbose && cat(verbose, "Resetting working directory: ", opwd);
    setwd(opwd);
  });


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup gatk arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The output directory should always to be the current directory
  opts <- c();

  if (!is.null(analysisType)) opts <- c(opts, "--analysis_type", shQuote(analysisType));
  if (!is.null(inputFile)) opts <- c(opts, "--input_file", shQuote(inputFile));
  if (!is.null(referenceSequence)) opts <- c(opts, "--reference_sequence", shQuote(referenceSequence));
  if (!is.null(interval)) opts <- c(opts, "--interval", shQuote(interval));
  opts <- c(opts, ...);
  str(opts);
  opts <- as.list(opts);
  str(opts);
  opts$stdout <- TRUE;
  opts$stderr <- TRUE;

  # Assert no duplicated options
  names <- names(opts);
  names <- names[nchar(names) > 0L];
  dups <- names[duplicated(names)];
  if (length(dups) > 0L) {
    throw("Duplicated options detected: ", paste(sQuote(dups), collapse=", "));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call GATK executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling systemGATK()");
  args <- opts;
  verbose && cat(verbose, "Arguments:");
  verbose && print(verbose, args);
  args$verbose <- less(verbose, 10);
  res <- do.call(systemGATK, args=args);
  status <- attr(res, "status"); if (is.null(status)) status <- 0L;
  verbose && cat(verbose, "Results:");
  verbose && str(verbose, res);
  verbose && cat(verbose, "Status:");
  verbose && str(verbose, status);
  verbose && exit(verbose);

  # Successful?
  if (status == 0L) {
    # If we get this far, assume it was all successful.
    # Allow the temporary output path to be renamed to the
    # intended output path instead of the "error" one.
    outPathFinal <- outPathOrg;
  }

  verbose && exit(verbose);

  res
}) # gatk()


############################################################################
# HISTORY:
# 2014-03-14 [HB]
# o Created from tophat.R.
############################################################################


