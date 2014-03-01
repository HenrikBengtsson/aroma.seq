###########################################################################/**
# @RdocDefault systemFastQC
#
# @title "Calls the FastQC executable"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments specifying FastQC command line switches.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemFastQC", "default", function(..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments '...':
  args <- list(...);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling FastQC executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  perl <- findPerl(verbose=less(verbose, 50));
  verbose && cat(verbose, "Perl executable: ", perl);

  fastq <- findFastQC(verbose=less(verbose, 50));
  verbose && cat(verbose, "FastQC executable: ", fastq);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call FastQC
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # WORKAROUND: 'fastqc' returns before the standard output is available
  # to R (at least on Windows). This results in an empty string to R.
  # By piping to a temporary output file instead, we can pool that file
  # for results.
  resfile <- tempfile();
  on.exit(file.remove(resfile), add=TRUE);
  args <- c(sprintf('"%s"', fastq), ...);
  verbose && cat(verbose, "Perl call: ", paste(args, collapse=" "));
  res <- system2(perl, args=args, stdout=resfile, stderr=resfile);
  verbose && cat(verbose, "Result code: ", res);

  # Pool file every 0.1 seconds for 5 seconds.
  for (kk in 1:50) {
    bfr <- tryCatch({
      suppressWarnings({
        readLines(resfile);
      })
    }, error = function(ex) {""})
    if (length(bfr) > 1L) break;
    if (length(bfr) == 1L && nchar(bfr[1L]) > 0L) break;
    Sys.sleep(0.1);
  }

  attr(bfr, "returnCode") <- res;

  verbose && exit(verbose);

  bfr;
}) # systemFastQC()


############################################################################
# HISTORY:
# 2014-02-28 [HB]
# o Created from systemPicard().
############################################################################
