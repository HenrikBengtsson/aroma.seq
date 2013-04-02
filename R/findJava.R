###########################################################################/**
# @RdocFunction findJava
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
findJava <- function(mustExists=TRUE, ..., verbose=FALSE) {
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

  # Validate by retrieving 'version' attribute.
  if (isFile(pathname)) {
    verbose && enter(verbose, "Retrieving version");
    res <- system2(pathname, "-version", stdout=TRUE, stderr=TRUE);
    ver <- grep("java version", res, value=TRUE);
    ver <- gsub("java version", "", ver);
    ver <- gsub("[ \"]", "", ver);
    ver <- gsub("_", "-", ver);
    ver <- package_version(ver);
    attr(pathname, "version") <- ver;
##    attr(pathname, "raw_version") <- res;
    verbose && exit(verbose);

    .findCache(name=command, path=pathname);
  } else if (mustExists) {
    throw(sprintf("Failed to located Java (executable '%s').", command));
  }

  verbose && exit(verbose);

  pathname;
} # findJava()


############################################################################
# HISTORY:
# 2013-04-01
# o Now findJava() sets attribute 'version', iff possible.
# 2012-09-28
# o Created.
############################################################################
