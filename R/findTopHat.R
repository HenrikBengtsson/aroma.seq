###########################################################################/**
# @RdocDefault findTopHat
# @alias findTopHat1
# @alias findTopHat1.default
# @alias findTopHat2
# @alias findTopHat2.default
#
# @title "Locates the TopHat executable of a certain version"
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
#   \item{command}{A @character string specifying the name of the executable.}
#   \item{version}{If non-@NULL, specifies which version of TopHat.}
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
setMethodS3("findTopHat", "default", function(mustExists=TRUE, command="tophat", version, ..., verbose=FALSE) {
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
}) # findTopHat()

setMethodS3("findTopHat1", "default", function(...) {
  findTopHat(..., command="tophat", version="1");
}) # findTopHat1()

setMethodS3("findTopHat2", "default", function(...) {
  res <- tryCatch({
    findTopHat(..., command="tophat", version="2");
  }, error = function(ex) { NULL });
  if (is.null(res)) {
    res <- findTopHat(..., command="tophat2", version="2");
  }
  res;
}) # findTopHat2()


############################################################################
# HISTORY:
# 2013-04-01
# o CLEANUP: Now findTopHat1() and  findTopHat2() utilizes findTopHat().
# o Now findTopHat() cached the results.
# o Renamed from findTopHatv() to findTopHat().
# 2013-01-24
# o Created.
############################################################################
