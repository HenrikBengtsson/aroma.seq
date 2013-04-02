###########################################################################/**
# @RdocDefault systemPicard
#
# @title "Calls the Picard executable"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{command}{A @character string specifying the Picard command
#     (the name of the jar file without the *.jar extension).}
#   \item{...}{Additional arguments specifying Picard command line switches.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{\dontrun{
#   pathnameBAM <- "bwaData/LambdaVirusExample,bwa,is/Generic/reads_1.bam"
#   res <- systemPicard("ValidateSamFile", INPUT=pathnameBAM, stderr=FALSE)
#   ## ERROR: Read groups is empty
#   print(res)
#   ## [1] 1
#
#   res <- systemPicard("ValidateSamFile", INPUT=pathnameBAM,
#                       IGNORE="MISSING_READ_GROUP", stderr=FALSE)
#   ## "No errors found"
#   print(res)
#   ## [1] 0
# }}
#
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemPicard", "default", function(command, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'command':
  command <- Arguments$getCharacter(command);

  # Arguments '...':
  args <- list(...);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling Picard executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- findPicard(verbose=less(verbose, 50));
  verbose && cat(verbose, "Picard directory: ", path);

  verbose && cat(verbose, "Picard command: ", command);
  filename <- sprintf("%s.jar", command);

  pathname <- file.path(path, filename);
  verbose && cat(verbose, "Pathname to Picard jar file: ", pathname);
  pathname <- Arguments$getReadablePathname(pathname);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call Picard java jar
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  .fmtArg <- c(".*"="%s=%s");
  res <- systemJavaJar(pathname, ..., .fmtArg=.fmtArg, verbose=less(verbose, 5));

  verbose && exit(verbose);

  res;
}) # systemPicard()


############################################################################
# HISTORY:
# 2012-09-28
# o Now utilizing systemJavaJar().
# 2012-09-27
# o Created.
############################################################################
