###########################################################################/**
# @RdocDefault findSamtools
#
# @title "Locates the samtools executable"
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
# @author
#*/###########################################################################
setMethodS3("findSamtools", "default", function(mustExists=TRUE, ..., verbose=FALSE) {
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

  verbose && enter(verbose, "Locating samtools software");

  command <- "samtools";
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
    verbose && enter(verbose, "Retrieving version");
    suppressWarnings({
      res <- system2(pathname, stdout=TRUE, stderr=TRUE);
    });
    ver <- grep("Version:", res, value=TRUE);
    ver <- gsub("Version:[ ]*([0-9.-]+).*", "\\1", ver);
    ver <- gsub("_", "-", ver);
    ver <- package_version(ver);
    attr(pathname, "version") <- ver;
    verbose && exit(verbose);
    .findCache(name=command, path=pathname);
  } else if (mustExists) {
    throw(sprintf("Failed to located samtools (executable '%s').", command));
  }

  verbose && exit(verbose);

  pathname;
}) # findSamtools()


############################################################################
# HISTORY:
# 2012-09-25
# o Created.
############################################################################
