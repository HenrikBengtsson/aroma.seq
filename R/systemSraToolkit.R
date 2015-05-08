###########################################################################/**
# @RdocDefault systemSraToolkit
#
# @title "Calls an SRA Toolkit executable"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{command}{The SRA Toolkit executable to call.}
#   \item{...}{Command line switches.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{\dontrun{
#   pathnameSRA <- "sraData/DataSetA/Homo_sapiens/reads.sra"
#   res <- systemSraToolkit("fastq-dump", pathnameSRA, ..., stderr=FALSE)
# }}
#
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemSraToolkit", "default", function(command, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'command':
  command <- Arguments$getFilename(command)

  # Arguments '...':
  args <- list(...);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling SRA Toolkit executable");
  verbose && cat(verbose, "Command: ", command);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- findSraToolkit(command, verbose=less(verbose, 50));
  verbose && cat(verbose, "SRA Toolkit executable: ", pathname);
  pathname <- Arguments$getReadablePathname(pathname);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call SRA Toolkit executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- system2(pathname, ...);

  verbose && exit(verbose);

  res;
}) # systemSraToolkit()


############################################################################
# HISTORY:
# 2014-09-29
# o Created.
############################################################################
