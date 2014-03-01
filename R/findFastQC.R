###########################################################################/**
# @RdocFunction findFastQC
# \alias{FASTQC_HOME}
#
# @title "Locates the FastQC executable"
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
#  The FastQC tools directory is searched for as follows:
#  \enumerate{
#   \item \code{Sys.getenv("FASTQC_HOME")}
#  }
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
findFastQC <- function(mustExist=TRUE, ..., verbose=FALSE) {
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

  verbose && enter(verbose, "Locating FastQC software");

  command <- "fastqc";
  verbose && cat(verbose, "Command: ", command);

  # Check for cached results
  res <- .findCache(name=command);
  if (!is.null(res)) {
    path <- res$path;
    verbose && cat(verbose, "Found cached result.");
    verbose && exit(verbose);
    return(path);
  }

  path <- Sys.getenv("FASTQC_HOME");
  verbose && printf(verbose, "System variable 'FASTQC_HOME': '%s'\n", path);
  if (path == "") path <- NULL;

  if (!is.null(path) && isDirectory(path)) {
    path <- Arguments$getReadablePath(path, mustWork=FALSE);
    verbose && cat(verbose, "Located directory: ", path);
    # Validating
    pathname <- file.path(path, "fastqc");
  } else {
    pathname <- Sys.which("fastqc");
  }

  # Done?
  if (!isFile(pathname)) {
    if (mustExist) {
      throw("Failed to located FastQC tools");
    } else {
      return(NULL);
    }
  }

  # Infer version
  verbose && enter(verbose, "Retrieving version");
  perl <- findPerl();

  # WORKAROUND: 'fastqc' returns before the standard output is available
  # to R (at least on Windows). This results in an empty string to R.
  # By piping to a temporary output file instead, we can pool that file
  # for results.
  resfile <- tempfile();
  on.exit(file.remove(resfile), add=TRUE);
  cmd <- sprintf('"%s" --version', pathname);
  res <- system2(perl, args=cmd, stdout=resfile);
  verbose && cat(verbose, "Result code: ", res);
  if (res != 0L) {
    throw("Failed to run 'fastqc --version'.");
  }

  # Pool file every 0.1 seconds for 5 seconds.
  for (kk in 1:50) {
    bfr <- tryCatch({
      suppressWarnings({
        readLines(resfile);
      })
    }, error = function(ex) {""})
    if (length(bfr) > 1L) break;
    if (length(bfr) == 1L && nchar(bfr[1L]) > 0L) break;
    Sys.sleep(0.1);
  }

  if (length(bfr) == 0L) {
    warning("Failed to infer FastQC version.");
  } else {
    bfr <- bfr[nchar(bfr) > 0L];
    ver <- bfr[1L];
    ver <- gsub("FastQC v([0-9.-_]+).*", "\\1", ver);
    # Try to coerce
    tryCatch({
      ver <- package_version(ver);
    }, error = function(ex) {})
  }
  attr(pathname, "version") <- ver;
  verbose && exit(verbose);

  .findCache(name=command, path=pathname);

  verbose && exit(verbose);

  pathname;
} # findFastQC()


############################################################################
# HISTORY:
# 2014-02-28
# o Created from findPicard().
############################################################################
