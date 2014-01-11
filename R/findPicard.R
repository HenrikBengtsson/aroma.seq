###########################################################################/**
# @RdocFunction findPicard
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
#   \item{mustExist}{If @TRUE, an exception is thrown if the executable
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
# @author "HB"
#
# @keyword internal
#*/###########################################################################
findPicard <- function(mustExist=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'mustExist':
  mustExist <- Arguments$getLogical(mustExist);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating Picard software");

  command <- "picard";
  verbose && cat(verbose, "Command: ", command);

  # Check for cached results
  res <- .findCache(name=command);
  if (!is.null(res)) {
    path <- res$path;
    verbose && cat(verbose, "Found cached result.");
    verbose && exit(verbose);
    return(path);
  }

  path <- Sys.getenv("PICARD_HOME");
  verbose && printf(verbose, "System variable 'PICARD_HOME': '%s'\n", path);
  if (path == "") path <- NULL;

  if (!is.null(path) && isDirectory(path)) {
    path <- Arguments$getReadablePath(path, mustWork=FALSE);
    verbose && cat(verbose, "Located directory: ", path);

    # Validating
    files <- list.files(path=path, pattern="[.]jar$");
    if (length(files) == 0L) {
      throw("The located Picard directory contains no *.jar files: ", path);
    }

    # Validate by retrieving 'version' attribute.
    pathnameT <- file.path(path, "ViewSam.jar");
    verbose && enter(verbose, "Retrieving version");
    res <- systemJavaJar(pathnameT, "--version", stdout=TRUE, stderr=TRUE);
    ver <- res[1L];
    ver <- gsub("([0-9.-_]+).*", "\\1", ver);
    # Try to coerce
    tryCatch({
      ver <- package_version(ver);
    }, error = function(ex) {})
    attr(path, "version") <- ver;
##    attr(pathname, "raw_version") <- res;
    verbose && exit(verbose);

    .findCache(name=command, path=path);
  } else if (mustExist) {
    throw(sprintf("Failed to located Picard tools"));
  }

  verbose && exit(verbose);

  path;
} # findPicard()


############################################################################
# HISTORY:
# 2013-04-01
# o Now findPicard() sets attribute 'version', iff possible.
# 2012-10-01
# o BUG FIX: findPicard(mustExists=FALSE) would throw an error if Picard
#   was not found.
# 2012-09-27
# o Created.
############################################################################
