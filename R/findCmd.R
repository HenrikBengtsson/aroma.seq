###########################################################################/**
# @RdocDefault findCmd
#
# @title "Locates the executable given by 'cmd'"
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
#  The executable is searched for as follows:
#  \enumerate{
#   \item \code{Sys.which("tophat")}
#  }
# }
#
# @author
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

  pathname <- Sys.which(command);
  if (identical(pathname, "")) pathname <- NULL;
  if (!isFile(pathname)) pathname <- NULL;
  verbose && cat(verbose, "Located pathname: ", pathname);
  if (mustExists && !isFile(pathname)) {
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
