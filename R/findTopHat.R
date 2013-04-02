###########################################################################/**
# @RdocFunction findTopHat
# @alias findTopHat1
# @alias findTopHat2
#
# @title "Locates an external executable"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{mustExists}{If @TRUE, an exception is thrown if the executable
#      could not be located.}
#   \item{command}{A @character string specifying the name of the
#      executable to locate.}
#   \item{version}{(optional) If non-@NULL, specifies which version of the
#      executable to trieve.}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the pathname (or the path) of the external executable.
#   If not found, @NULL is returned, unless if \code{mustExists=TRUE}
#   in case an error is thrown.
# }
#
# \details{
#  The TopHat executable is searched for as follows:
#  \enumerate{
#   \item \code{Sys.which(command)}
#  }
# }
#
# @author
#*/###########################################################################
findTopHat <- function(mustExists=TRUE, command="tophat", version, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'mustExists':
  mustExists <- Arguments$getLogical(mustExists);

  # Argument 'command':
  command <- Arguments$getCharacter(command);

  # Argument 'version':
  version <- Arguments$getCharacter(version);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating TopHat software");

  verbose && cat(verbose, "Command: ", command);
  verbose && cat(verbose, "Version: ", version);

  # Check for cached results
  res <- .findCache(name=command, version=version);
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
    # Check version
    res <- system2(pathname, args="--version", stdout=TRUE);
    ver <- sub("^[^0-9]+", "", res);
    ver <- package_version(ver);
    # Make sure that TopHat v1 was retrieved
    if (ver < version || ver > version) {
      throw(sprintf("Failed to located TopHat v%s. Found TopHat v%s", version, ver));
    }
    .findCache(name=command, version=version, path=pathname);
  } else if (mustExists) {
    throw(sprintf("Failed to located TopHat (executable '%s').", command));
}

  verbose && exit(verbose);

  pathname;
} # findTopHat()


findTopHat1 <- function(...) {
  findTopHat(..., command="tophat", version="1");
} # findTopHat1()


findTopHat2 <- function(...) {
  res <- tryCatch({
    findTopHat(..., command="tophat", version="2");
  }, error = function(ex) { NULL });
  if (is.null(res)) {
    res <- findTopHat(..., command="tophat2", version="2");
  }
  res;
} # findTopHat2()


############################################################################
# HISTORY:
# 2013-04-01
# o CLEANUP: Now findTopHat1() and  findTopHat2() utilizes findTopHat().
# o Now findTopHat() cached the results.
# o Renamed from findTopHatv() to findTopHat().
# 2013-01-24
# o Created.
############################################################################
