###########################################################################/**
# @RdocDefault findCmd
#
# @title "Locates the executable given by 'command'"
#
# \description{
#  @get "title" on the current system.
# }
#
# @synopsis
#
# \arguments{
#   \item{command}{Name of executable to find}
#   \item{mustExists}{If @TRUE, an exception is thrown if the executable
#      could not be located.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \details{
#  The executable is searched for as follows:
#  \enumerate{
#   \item \code{Sys.which(command)}
#  }
#  NB: This method does NOT do any version checking!
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("findCmd", "default", function(command, mustExists=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'command':
  command <- Arguments$getCharacter(command);
  if (nchar(command) < 1) {
    throw(sprintf("Failed to locate executable:  command argument is NULL."));
  }

  # Argument 'mustExists':
  mustExists <- Arguments$getLogical(mustExists);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating software");
  verbose && cat(verbose, "Command: ", command);

  # Check for cached results
  res <- .findCache(name=command);
  if (!is.null(res)) {
    pathname <- res$path;
    verbose && cat(verbose, "Found cached result.");
    verbose && exit(verbose);
    return(pathname);
  }

  pathname <- Sys.which(command);
  if (identical(pathname, "")) pathname <- NULL;
  if (!isFile(pathname)) pathname <- NULL;
  verbose && cat(verbose, "Located pathname: ", pathname);

  if (isFile(pathname)) {
    .findCache(name=command, path=pathname);
  } else if (mustExists) {
    throw(sprintf("Failed to locate (executable '%s').", command));
  }

  verbose && exit(verbose);

  pathname;
}) # findCmd()


############################################################################
# HISTORY:
# 2013-01-29
# o Created TAT
############################################################################
