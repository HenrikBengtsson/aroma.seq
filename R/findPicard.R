###########################################################################/**
# @RdocDefault findPicard
# \alias{PICARD_HOME}
#
# @title "Locates the Picard executable"
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
#  The Picard tools directory is searched for as follows:
#  \enumerate{
#   \item \code{Sys.getenv("PICARD_HOME")}
#  }
# }
#
# @author
#*/###########################################################################
setMethodS3("findPicard", "default", function(mustExists=TRUE, ..., verbose=FALSE) {
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

  verbose && enter(verbose, "Locating Picard software");

  command <- "picard";
  verbose && cat(verbose, "Command: ", command);

  path <- Sys.getenv("PICARD_HOME");
  verbose && printf(verbose, "System variable 'PICARD_HOME': '%s'\n", path);
  if (path == "") path <- NULL;

  if (!is.null(path) && isDirectory(path)) {
    path <- Arguments$getReadablePath(path);
    verbose && cat(verbose, "Located directory: ", path);

    # Validating
    files <- list.files(path=path, pattern="[.]jar$");
    if (length(files) == 0L) {
      throw("The located Picard directory contains no *.jar files: ", path);
    }
  }
 
  if (mustExists && (is.null(path) || !isDirectory(path))) {
    throw(sprintf("Failed to located Picard tools"));
  }

  verbose && exit(verbose);

  path;
}) # findPicard()


############################################################################
# HISTORY:
# 2012-10-01
# o BUG FIX: findPicard(mustExists=FALSE) would throw an error if Picard
#   was not found.
# 2012-09-27
# o Created.
############################################################################
