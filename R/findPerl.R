###########################################################################/**
# @RdocFunction findPerl
#
# @title "Locates the Perl executable"
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
#  The Perl executable is searched for as follows:
#  \enumerate{
#   \item \code{Sys.which("perl")}
#  }
# }
#
# @author
#*/###########################################################################
findPerl <- function(mustExists=TRUE, ..., verbose=FALSE) {
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

  verbose && enter(verbose, "Locating Perl software");

  command <- "perl";
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
    ver <- grep("This is perl", res, value=TRUE);
    ver <- gsub(".*[(]v([0-9.]+)[)].*", "\\1", ver);
    ver <- gsub("[ \"]", "", ver);
    ver <- gsub("_", "-", ver);
    ver <- package_version(ver);
    attr(pathname, "version") <- ver;
##    attr(pathname, "raw_version") <- res;
    verbose && exit(verbose);

    .findCache(name=command, path=pathname);
  } else if (mustExists) {
    throw(sprintf("Failed to located Perl (executable '%s').", command));
  }

  verbose && exit(verbose);

  pathname;
} # findPerl()


############################################################################
# HISTORY:
# 2013-04-01
# o Now findPerl() sets attribute 'version', iff possible.
# 2012-09-28
# o Created.
############################################################################
