###########################################################################/**
# @RdocDefault findJava
#
# @title "Locates the Java executable"
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
#  The Java executable is searched for as follows:
#  \enumerate{
#   \item \code{Sys.which("java")}
#  }
# }
#
# @author
#*/###########################################################################
setMethodS3("findJava", "default", function(mustExists=TRUE, ..., verbose=FALSE) {
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

  verbose && enter(verbose, "Locating Java software");

  command <- "java";
  verbose && cat(verbose, "Command: ", command);

  pathname <- Sys.which(command);
  if (identical(pathname, "")) pathname <- NULL;
  if (!isFile(pathname)) pathname <- NULL;


  verbose && cat(verbose, "Located pathname: ", pathname);

  if (mustExists && !isFile(pathname)) {
    throw(sprintf("Failed to located Java (executable '%s').", command));
  }

  verbose && exit(verbose);

  pathname;
}) # findJava()


############################################################################
# HISTORY:
# 2012-09-28
# o Created.
############################################################################
