###########################################################################/**
# @RdocDefault systemGATK
#
# @title "Calls the GATK executable"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments specifying GATK command line switches.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{\dontrun{
#   pathnameBAM <- "bwaData/LambdaVirusExample,bwa,is/Generic/reads_1.bam"
#   res <- systemGATK(T="CountReads", ..., stderr=FALSE)
# }}
#
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemGATK", "default", function(..., verbose=FALSE) {
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

  verbose && enter(verbose, "Calling GATK executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- findGATK(verbose=less(verbose, 50));
  verbose && cat(verbose, "GATK jar file: ", pathname);
  pathname <- Arguments$getReadablePathname(pathname);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call GATK java jar
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  .fmtArg <- c("--.*"="%s %s", ".*"="-%s %s");
  res <- systemJavaJar(pathname, ..., .fmtArg=.fmtArg, verbose=less(verbose, 5));

  verbose && exit(verbose);

  res;
}) # systemGATK()


############################################################################
# HISTORY:
# 2012-09-28
# o Created from systemPicard.R.
############################################################################
