###########################################################################/**
# @RdocFunction findBowtie2
#
# @title "Locates one of the bowtie2 executable"
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
#   \item{command}{A @character string specifying the executable to locate.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \details{
#  The Bowtie2 executable is searched for as follows:
#  \enumerate{
#   \item \code{Sys.which(command)}
#  }
# }
#
# @author
#*/###########################################################################
findBowtie2 <- function(mustExists=TRUE, ..., command=c("bowtie2", "bowtie2-align", "bowtie2-build", "bowtie2-inspect"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'mustExists':
  mustExists <- Arguments$getLogical(mustExists);

  # Argument 'command':
  command <- match.arg(command);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating Bowtie2 software");
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
    res <- system2(pathname, args="--version", stdout=TRUE, stderr=TRUE);
    ver <- grep("version", res, value=TRUE)[1L];
    ver <- gsub(".*version[ ]*([0-9.-_]+.*)", "\\1", ver);
    ver <- gsub("_", "-", ver);
    # Try to coerce
    tryCatch({
      ver <- package_version(ver);
    }, error = function(ex) {})
    attr(pathname, "version") <- ver;
##    attr(pathname, "raw_version") <- res;
    verbose && exit(verbose);

    .findCache(name=command, path=pathname);
  } else if (mustExists) {
    throw(sprintf("Failed to located Bowtie2 (executable '%s').", command));
  }

  verbose && exit(verbose);

  pathname;
} # findBowtie2()

############################################################################
# HISTORY:
# 2013-04-01
# o BUG FIX: findBowtie2() was incorrectly hardwired to 'bowtie2-align'.
# o Now findBowtie2() sets attribute 'version', iff possible.
# 2012-09-27
# o Now looking for bowtie2-align instead of bowtie2, because the latter
#   is a perl script calling the former.
# 2012-09-24
# o Created.
############################################################################
