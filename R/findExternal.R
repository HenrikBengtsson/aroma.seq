###########################################################################/**
# @RdocFunction findExternal
# @alias getExternalHome
# @alias findJava
# @alias findPerl
# @alias findPython
# @alias findBowtie2
# @alias findBWA
# @alias findHTSeq
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
#   # The generic internal function used
#   @usage findExternal
#
#   # Programming environments
#   @usage findJava
#   @usage findPerl
#   @usage findPython
#
#   # Samtools
#   @usage findSamtools
#
#   # HTSeq
#   @usage findHTSeq
#
#   # BWA
#   @usage findBWA
#
#   # Bowtie and TopHat
#   @usage findBowtie2
#   @usage findTopHat
#   @usage findTopHat1
#   @usage findTopHat2
#
#   # HTSeq
#   @usage findHTSeq
#
#   # fastq-dump
#   @usage findFastqDump
# }
#
# \arguments{
#   \item{mustExist}{If @TRUE, an exception is thrown if the executable
#      could not be located.}
#   \item{command}{A @character string specifying the name of the
#      executable to locate.}
#   \item{version}{If non-@NULL, specifies which version of the
#      executable to retrieve.}
#   \item{versionPattern}{(A named @character string specifying the
#      @see "base::gsub" regular expression to extraction the version
#      where there name is the command-line option specifying how
#      to call the external for retrieving the version output.}
#   \item{expectedStatus}{An @integer @vector of expected status codes
#      returned when querying the executable for the version.}
#   \item{force}{If @TRUE, cached results are ignored, otherwise not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{...}{Additional arguments passed to @see "findExternal", or ignored.}
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
findExternal <- function(mustExist=TRUE, command, path=NULL, version=NULL, versionPattern=NULL, expectedStatus=c(0L, 1L), force=FALSE, verbose=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'mustExist':
  mustExist <- Arguments$getLogical(mustExist);

  # Argument 'command':
  command <- Arguments$getCharacter(command);

  # Argument 'path':
  if (!is.null(path)) {
    path <- Arguments$getReadablePath(path);
  }

  # Argument 'version':
  if (!is.null(version)) {
    version <- Arguments$getCharacters(version);
    if (length(version) == 1L) {
      version <- rep(version, length.out=2L);
    }
  }

  # Argument 'versionPattern':
  if (!is.null(versionPattern)) {
    name <- names(versionPattern);
    versionPattern <- Arguments$getRegularExpression(versionPattern);
    names(versionPattern) <- name;
  } else if (!is.null(version)) {
    throw("Argument 'versionPattern' must be specified if 'version' is: ", version);
  }

  # Argument 'expectedStatus':
  expectedStatus <- Arguments$getIntegers(expectedStatus);

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
    verbose && printf(verbose, "Requested version range: [%s,%s)\n", version[1L], version[2L]);
    verbose && cat(verbose, "Version option: ", verOpt);
    verbose && cat(verbose, "Version pattern: ", versionPattern);
  }

  # Check for cached results
  if (!force) {
    res <- .findCache(name=command, version=version, versionPattern=versionPattern);
    if (!is.null(res)) {
      pathname <- res$path;
      if (!is.null(pathname)) {
        verbose && cat(verbose, "Found cached result.");
        verbose && exit(verbose);
        return(pathname);
      }
    }
  }


  pathname <- NULL;

  # (a) Search in predefined 'path'?
  if (!is.null(path)) {
    # Search executables in the specified 'path'...
    filenames <- NULL;
    if (.Platform$OS.type == "windows") {
      # Let *.bat and *.exe override the ones without extensions
      filenames <- c(filenames, sprintf("%s.bat", command));
      filenames <- c(filenames, sprintf("%s.exe", command));
    }
    # Let *.sh override the ones without extensions
    filenames <- c(filenames, sprintf("%s.sh", command));
    filenames <- c(filenames, command);
    for (filename in filenames) {
      pathname <- file.path(path, filename);
      if (isFile(pathname)) break;
      pathname <- NULL;
    }
  }


  # (b) Search system 'PATH'?
  if (is.null(pathname)) {
    pathname <- Sys.which(command);
    if (identical(pathname, "")) pathname <- NULL;
  }
  verbose && cat(verbose, "Located pathname: ", pathname);


  # (c) Is the executable required?
  if (mustExist && !isFile(pathname)) {
    throw(sprintf("Failed to locate external executable '%s'", command));
  }


  # (d) Query version?
  if (!is.null(versionPattern) && isFile(pathname)) {
    verbose && enter(verbose, "Retrieving version");

    # Request version output from software
    suppressWarnings({
      res <- system2(pathname, args=verOpt, stdout=TRUE, stderr=TRUE);
    });

    # Status code
    status <- attr(res, "status");
    verbose && cat(verbose, "Return status: ", status);

    # Validate return status code
    if (length(status) > 0L && length(expectedStatus) > 0L) {
      if (!is.element(status, expectedStatus)) {
        throw(sprintf("Unexpected return status code when calling %s: %d != (%s)",
                      sQuote(pathname), status, paste(expectedStatus, collapse=", ")));
      }
    }

    # Parse
    resT <- paste(res, collapse=" ");  # Search across newlines
    ver <- grep(versionPattern, resT, value=TRUE);
    if (length(ver) > 0L) {
      ver <- ver[1L];
      verbose && printf(verbose, "Version (output): '%s'\n", ver);
      ver <- gsub(versionPattern, "\\1", ver);
      verbose && printf(verbose, "Version (string): '%s'\n", ver);
      # Drop trailing periods and more
      ver <- gsub("[.]$", "", ver);
      ver <- trim(ver);
      verbose && printf(verbose, "Version (trimmed): '%s'\n", ver);
      # Try to coerce
      tryCatch({
        ver <- gsub("_", "-", ver);
        ver <- package_version(ver);
        verbose && printf(verbose, "Version (parsed): '%s'\n", ver);
      }, error = function(ex) {});
    } else {
      msg <- sprintf("Failed to identify 'version' using regular expression '%s': %s", versionPattern, paste(res, collapse="\\n"));
      if (!is.null(version)) {
        throw(msg);
      } else {
        warning(msg);
      }
      ver <- NULL;
    }

    verbose && exit(verbose);


    # Record the version
    attr(pathname, "version") <- ver;

    if (!is.null(version)) {
      verbose && enter(verbose, "Validated version");
      verbose && cat(verbose, "Available version: ", ver);
      verbose && printf(verbose, "Requested version range: [%s,%s)\n", version[1L], version[2L]);
      if (ver < version[1L] || ver >= version[2L]) {
        pathname <- NULL;
        if (mustExist) {
          throw(sprintf("Failed to locate external executable '%s' with version in [%s,%s): %s", command, version[1L], version[2L], ver));
        }
      }
      verbose && exit(verbose);
    }
  } # if (!is.null(versionPattern) && isFile(pathname))


  # Save cache
  .findCache(name=command, version=version, versionPattern=versionPattern, path=pathname);

  verbose && exit(verbose);

  pathname;
} # findExternal()


getExternalHome <- function(name, mustExist=TRUE, ...) {
  path <- Sys.getenv(name);
  if (path == "") path <- NULL;
  if (!is.null(path) && mustExist) {
    path <- Arguments$getReadablePath(path, mustExist=FALSE);
    if (!isDirectory(path)) {
      throw(sprintf("System environment variable %s specifies a non-existing directory: %s", sQuote(name), path));
    }
  }
  path;
} # getExternalHome()


############################################################################
# HISTORY:
# 2014-08-07
# o Added argument 'path' to findExternal().
# o Added getExternalHome().
# 2014-07-24
# o ROBUSTNESS: Added 'versionPattern' to the set of cache keys.
# 2014-03-09
# o Now findExternal() drops trailing periods and trims too.
# 2013-04-02
# o ROBUSTNESS: Now findExternal(..., version=NULL) only gives a warning
#   if it fails to infer the version from the software's version output.
# 2013-04-01
# o Renamed argument 'mustExists' to 'mustExist'.
# o Created from findTopHat.R.
############################################################################
