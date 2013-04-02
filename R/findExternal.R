###########################################################################/**
# @RdocFunction findExternal
# @alias findJava
# @alias findPerl
# @alias findBowtie2
# @alias findBWA
# @alias findSamtools
# @alias findTopHat
# @alias findTopHat1
# @alias findTopHat2
#
# @title "Locates an external executable"
#
# \description{
#  @get "title".
# }
#
# \usage{
#   findJava(@eval "t<-formals(findJava);paste(gsub('=$','',paste(names(t),t,sep='=')),collapse=', ')")
#   findPerl(@eval "t<-formals(findPerl);paste(gsub('=$','',paste(names(t),t,sep='=')),collapse=', ')")
#   findBowtie2(@eval "t<-formals(findBowtie2);paste(gsub('=$','',paste(names(t),t,sep='=')),collapse=', ')")
#   findBWA(@eval "t<-formals(findBWA);paste(gsub('=$','',paste(names(t),t,sep='=')),collapse=', ')")
#   findSamtools(@eval "t<-formals(findSamtools);paste(gsub('=$','',paste(names(t),t,sep='=')),collapse=', ')")
#   findTopHat(@eval "t<-formals(findTopHat);paste(gsub('=$','',paste(names(t),t,sep='=')),collapse=', ')")
#   findTopHat1(@eval "t<-formals(findTopHat1);paste(gsub('=$','',paste(names(t),t,sep='=')),collapse=', ')")
#   findTopHat2(@eval "t<-formals(findTopHat2);paste(gsub('=$','',paste(names(t),t,sep='=')),collapse=', ')")
# }
#
# \arguments{
#   \item{mustExist}{If @TRUE, an exception is thrown if the executable
#      could not be located.}
#   \item{command}{A @character string specifying the name of the
#      executable to locate.}
#   \item{version}{If non-@NULL, specifies which version of the
#      executable to trieve.}
#   \item{versionPattern}{(A named @character string specifying the
#      @see "base::gsub" regular expression to extraction the version
#      where there name is the command-line option specifying how
#      to call the external for retrieving the version output.}
#   \item{force}{If @TRUE, cached results are ignored, otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{...}{Not used.}
# }
#
# \value{
#   Returns the pathname (or the path) of the external executable.
#   If not found, @NULL is returned, unless if \code{mustExist=TRUE}
#   in case an error is thrown.
#   If \code{versionPattern} is specified, then the inferred version
#   is returned as attribute 'version'.
# }
#
# \details{
#  The executable is searched using (in order):
#  \enumerate{
#   \item \code{Sys.which(command)}
#  }
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
findExternal <- function(mustExist=TRUE, command, version=NULL, versionPattern=NULL, force=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'mustExist':
  mustExist <- Arguments$getLogical(mustExist);

  # Argument 'command':
  command <- Arguments$getCharacter(command);

  # Argument 'version':
  if (!is.null(version)) {
    version <- Arguments$getCharacters(version);
    version <- rep(version, length.out=2L);
  }

  # Argument 'versionPattern':
  if (!is.null(versionPattern)) {
    name <- names(versionPattern);
    versionPattern <- Arguments$getRegularExpression(versionPattern);
    names(versionPattern) <- name;
  } else if (!is.null(version)) {
    throw("Argument 'versionPattern' must be specified if 'version' is: ", version);
  }

  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Locating external software");
  verOpt <- names(versionPattern);
  verbose && cat(verbose, "Command: ", command);
  if (!is.null(version)) {
    verbose && printf(verbose, "Requested version range: [%s,%s]\n", version[1L], version[2L]);
    verbose && cat(verbose, "Version option: ", verOpt);
    verbose && cat(verbose, "Version pattern: ", versionPattern);
  }

  # Check for cached results
  if (!force) {
    res <- .findCache(name=command, version=version);
    if (!is.null(res)) {
      pathname <- res$path;
      verbose && cat(verbose, "Found cached result.");
      verbose && exit(verbose);
      return(pathname);
    }
  }

  pathname <- Sys.which(command);
  if (identical(pathname, "")) pathname <- NULL;
  if (!isFile(pathname)) pathname <- NULL;
  verbose && cat(verbose, "Located pathname: ", pathname);

  if (isFile(pathname)) {
    if (!is.null(versionPattern)) {
      verbose && enter(verbose, "Retrieving version");
      suppressWarnings({
        res <- system2(pathname, args=verOpt, stdout=TRUE, stderr=TRUE);
      });
      ver <- grep(versionPattern, res, value=TRUE);
      if (length(ver) == 0L) {
        throw(sprintf("Failed to identify 'version' using regular expression '%s': %s", versionPattern, paste(res, collapse="\\n")));
      }
      verbose && printf(verbose, "Version (output): '%s'\n", ver);
      ver <- gsub(versionPattern, "\\1", ver);
      verbose && printf(verbose, "Version (string): '%s'\n", ver);
      tryCatch({
        ver <- gsub("_", "-", ver);
        ver <- package_version(ver);
        verbose && printf(verbose, "Version (parsed): '%s'\n", ver);
      }, error = function(ex) {});
      verbose && exit(verbose);

      # Record the version
      attr(pathname, "version") <- ver;

      if (!is.null(version)) {
        verbose && enter(verbose, "Validated version");
        verbose && cat(verbose, "Available version: ", ver);
        verbose && printf(verbose, "Requested version range: [%s,%s]\n", version[1L], version[2L]);
        if (ver < version[1L] || ver > version[2L]) {
          throw(sprintf("Failed to located '%s' with version in [%s,%s]: %s", command, version[1L], version[2L], ver));
        }
        verbose && exit(verbose);
      }
    }
    .findCache(name=command, version=version, path=pathname);
  } else if (mustExist) {
    throw(sprintf("Failed to located external executable '%s'", command));
  }

  verbose && exit(verbose);

  pathname;
} # findExternal()


############################################################################
# HISTORY:
# 2013-04-01
# o Renamed argument 'mustExists' to 'mustExist'.
# o Created from findTopHat.R.
############################################################################
