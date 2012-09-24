###########################################################################/**
# @RdocDefault findBWA
#
# @title "Locates the BWA executable"
#
# \description{
#  @get "title" on the current system.
# }
# 
# @synopsis
#
# \arguments{
#   \item{mustExists}{If @TRUE, an exception is thrown if the executable
#      could not be located.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author
#*/###########################################################################
setMethodS3("findBWA", "default", function(..., mustExists=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'mustExists':
  mustExists <- Arguments$getLogical(mustExists);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating BWA software");

  command <- "bwa";
  verbose && cat(verbose, "Command: ", command);

  pathname <- Sys.which(command);
  if (identical(pathname, "")) pathname <- NULL;
  if (!isFile(pathname)) pathname <- NULL;


  verbose && cat(verbose, "Located pathname: ", pathname);

  if (mustExists && !isFile(pathname)) {
    throw(sprintf("Failed to located BWA (executable '%s').", command));
  }

  verbose && exit(verbose);

  pathname;
}) # findBWA()


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################
