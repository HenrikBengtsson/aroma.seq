###########################################################################/**
# @RdocDefault fastqDump
#
# @title "Calls the fastq-dump executable"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to @see "systemFastqDump".}
#   \item{outPath}{Directory where output fastq files will reside at completion}
#   \item{pathnames}{Zero or more input pathnames.}
#   \item{split3}{If TRUE, add '--split-3' option to produce two files in case of paired-end data}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# \references{
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("fastqDump", "default", function(...,
                                             split3=TRUE,
                                             outPath="fastqData",
                                             pathnames=character(0L), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathnames':
  pathnames <- lapply(pathnames, FUN=Arguments$getReadablePathname);
  pathnames <- unlist(pathnames, use.names=TRUE);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running FastqDump");
  verbose && cat(verbose, "Input pathnames:");
  verbose && print(verbose, pathnames);

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate output atomically by writing to a temporary directory
  # that is renamed upon completion.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  outPathOrg <- outPath;
  outPath <- sprintf("%s.TMP", outPath);
  outPath <- Arguments$getWritablePath(outPath, mustNotExist=TRUE);
  verbose && cat(verbose, "Temporary output directory: ", outPath);
  # At the end, assume failure, unless successful (whereupon outPathFinal is changed to outPathOrg)
  outPathFinal <- sprintf("%s.ERROR", outPath);
  onExit({
    # removeDirectory(outPathOrg, recursive=FALSE, mustExist=FALSE);
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
  
  
  args <- list(...);
  if (length(pathnames) > 0L) {
    # Always quote pathnames (to handle spaces)
    pathnames <- sprintf('"%s"', pathnames);
    args <- c(list(pathnames), args);
  }
  if (!is.null(outPath)) {
    args <- c(args, sprintf('-O "%s"', outPath))
  }
  if (split3) {
    args <- c(args, "--split-3")    
  }

  browser()
  
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, args);
  args <- c(args,verbose=less(verbose, 10));
  
  browser()
  
  res <- do.call(systemFastqDump, args=args);
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
  
  browser()
  
  verbose && exit(verbose);

  res;
}) # fastqDump()


############################################################################
# HISTORY:
# 2014-02-28
# o Created from bwaSamse().
############################################################################
