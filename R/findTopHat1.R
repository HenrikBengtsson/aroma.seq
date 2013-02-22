###########################################################################/**
# @RdocDefault findTopHat1
#
# @title "Locates the TopHat executable and tests version"
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
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \details{
#  The TopHat executable is searched for as follows:
#  \enumerate{
#   \item \code{Sys.which("tophat")}
#  }
# }
#
# @author
#*/###########################################################################

setMethodS3("findTopHat1", "default", function(mustExists=TRUE, ..., verbose=FALSE) {
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

  verbose && enter(verbose, "Locating TopHat software");

  command <- "tophat";
  verbose && cat(verbose, "Command: ", command);

  pathname <- Sys.which(command);
  if (identical(pathname, "")) pathname <- NULL;
  if (!isFile(pathname)) pathname <- NULL;


  verbose && cat(verbose, "Located pathname: ", pathname);

  if (mustExists && !isFile(pathname)) {
    throw(sprintf("Failed to located TopHat (executable '%s').", command));
  }

  versionStr <- system2(pathname, args="--version", stdout=TRUE)
  versionStr <- sub("^[^0-9]+", "", versionStr)
  stopifnot(package_version(versionStr) >= 1 && package_version(versionStr) < 2)

  verbose && exit(verbose);

  pathname;
}) # findTopHat()


############################################################################
# HISTORY:
# 2013-01-24
# o Created.
############################################################################
