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
#   \item{analysisType}{(required) A @character string specifying type of
#     analysis.
#     In GATK, this corrensponds to engine argument
#     \code{--analysis_type} (or \code{-T}).}
#   \item{...}{(optional) named arguments}.
#   \item{pathnameI}{(optional) A @character @vector specifying one or more
#     (SAM or BAM) files.
#     In GATK, this corrensponds to one or more engine arguments
#     \code{--input_file} (or \code{-I}).}
#   \item{pathnameR}{(optional) A @character string specifying a reference file.
#     In GATK, this corrensponds to engine argument
#     \code{--reference_sequence} (or \code{-R}).}
#   \item{pathnameL}{(optional) A @character string specifying a ROD file.
#     In GATK, this corrensponds to engine argument
#     \code{--interval} (or \code{-L}).}
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
setMethodS3("gatk", "default", function(analysisType, ..., pathnameI=NULL, pathnameR=NULL, pathnameL=NULL, pathnameLog=NULL, outPath=NULL, overwrite=FALSE, verbose=FALSE) {
  # Make sure to evaluate registered onExit() statements
  on.exit(eval(onExit()));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'analysisType':
  analysisType <- Arguments$getCharacter(analysisType);

  # Argument 'pathnameI':
  if (length(pathnameI) > 0L) {
    pathnameI <- Arguments$getReadablePathnames(pathnameI, absolutePath=TRUE);
    assertNoDuplicated(pathnameI);
  }

  # Argument 'pathnameR':
  if (length(pathnameR) > 0L) {
    pathnameR <- Arguments$getReadablePathname(pathnameR, absolutePath=TRUE);
  }

  # Argument 'pathnameL':
  if (length(pathnameL) > 0L) {
    pathnameL <- Arguments$getReadablePathname(pathnameL, absolutePath=TRUE);
  }

  # Argument 'pathnameLog':
  if (length(pathnameLog) > 0L) {
    pathnameLog <- Arguments$getWritablePathname(pathnameLog, absolutePath=TRUE);
  }

  # Argument 'outPath':
  if (is.null(outPath)) {
    checksum <- getChecksum(list(analysisType=analysisType, ..., pathnameI=pathnameI, pathnameR=pathnameR, pathnameL=pathnameL));
    outPath <- file.path("gatkData", checksum);
    outPath <- Arguments$getWritablePath(outPath, mustNotExist=!overwrite);
    removeDirectory(outPath, recursive=overwrite, mustExist=TRUE);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
     pushState(verbose);
     on.exit(popState(verbose), add=TRUE);
  }

  verbose && enter(verbose, "Running gatk()");
  verbose && cat(verbose, "Analysis type: ", analysisType);
  verbose && cat(verbose, "Output directory: ", outPath);

  if (length(pathnameI) > 0L) {
    verbose && cat(verbose, "Input files:");
    verbose && print(verbose, sapply(pathnameI, FUN=getRelativePath));
  }

  if (length(pathnameR) > 0L) {
    verbose && cat(verbose, "Reference file:");
    verbose && print(verbose, getRelativePath(pathnameR));
  }

  if (length(pathnameL) > 0L) {
    verbose && cat(verbose, "Interval ROD file:");
    verbose && print(verbose, getRelativePath(pathnameL));
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate output atomically by writing to a temporary directory
  # that is renamed upon completion.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  outPathT <- sprintf("%s.TMP", outPath);
  outPathT <- Arguments$getWritablePath(outPathT, mustNotExist=TRUE);
  verbose && cat(verbose, "Temporary output directory: ", outPathT);
  # At the end, assume failure, unless successful.
  outPathFinal <- sprintf("%s.ERROR", outPath);
  onExit({
    # Rename temporary output directory to the final one
    file.rename(outPathT, outPathFinal);
  });


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # While running, let the output directory be the working directory.
  # Make sure the working directory is restored when exiting.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opwd <- setwd(outPathT);
  verbose && cat(verbose, "Original working directory: ", opwd);
  onExit({
    verbose && cat(verbose, "Resetting working directory: ", opwd);
    setwd(opwd);
  });


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup gatk arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The output directory should always to be the current directory
  opts <- c("--analysis_type", shQuote(analysisType));
  for (p in pathnameI) opts <- c(opts, "--input_file", shQuote(p));
  for (p in pathnameR) opts <- c(opts, "--reference_sequence", shQuote(p));
  for (p in pathnameL) opts <- c(opts, "--interval", shQuote(p));
  for (p in pathnameLog) opts <- c(opts, "--log_to_file", shQuote(p));
  opts <- c(opts, ...);
  opts <- as.list(opts);
  opts$stdout <- TRUE;
  opts$stderr <- TRUE;


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
    outPathFinal <- outPath;
  } else {
    # Any errors?
    pattern <- "^(#)+( )+ERROR MESSAGE:( )+";
    errors <- grep(pattern, res, value=TRUE);
    if (length(errors) > 0L) {
      errors <- gsub(pattern, "", errors);
      errors <- trim(errors);
      errors <- paste(errors, collapse="\n");
      throw(errors);
    }
  }

  verbose && exit(verbose);

  res
}) # gatk()


############################################################################
# HISTORY:
# 2014-04-14 [HB]
# o Added argument 'pathnameLog' to gatk().
# o ROBUSTNESS: Now gatk() outputs to a "unique" directory.
# 2014-04-13 [HB]
# o ROBUSTNESS: Now gatk() throws GATK errors as R errors.
# o Added optional arguments 'pathnameI', 'pathnameR' and 'pathnameL'
#   to gatk().
# o Made argument 'analysisType' mantadory for gatk().
# 2014-03-14 [HB]
# o Created from tophat.R.
############################################################################
