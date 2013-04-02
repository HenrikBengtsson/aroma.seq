###########################################################################/**
# @RdocDefault findCmdv
#
# @title "Locates the executable given by 'command'; tests version if possible"
#
# \description{
#  @get "title" on the current system.
# }
#
# @synopsis
#
# \arguments{
#   \item{command}{Name of executable to find}
#   \item{mustExist}{If @TRUE, an exception is thrown if the executable
#      could not be located.}
#   \item{version}{}
#   \item{...}{Not used.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \details{
#  The executable is searched for as follows:
#  \enumerate{
#   \item \code{Sys.which(command)}
#  }
#  The software version is obtained by trying to parse
#  \enumerate{
#   \item \code{'cmd --version'}
#  }
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("findCmdv", "default", function(command=NULL,
                                            version=NULL,
                                            mustExist=TRUE, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'command':
  command <- Arguments$getCharacter(command);
  if (nchar(command) < 1) {
    throw(sprintf("Failed to locate executable:  command argument is NULL."));
  }

  # Argument 'mustExist':
  mustExist <- Arguments$getLogical(mustExist);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  # Argument 'version'"
  if (is.null(version)) {
    throw("findCmdv:  version argument is null")
  }

  verbose && enter(verbose, "Locating software");
  verbose && cat(verbose, "Command: ", command);

  pathname <- Sys.which(command);
  if (identical(pathname, "")) pathname <- NULL;
  if (!isFile(pathname)) pathname <- NULL;
  verbose && cat(verbose, "Located pathname: ", pathname);
  if (mustExist && !isFile(pathname)) {
    throw(sprintf("Failed to locate (executable '%s').", command));
  }

  ## Primitive parsing to get version string
  versionStr <- system2(pathname, args="--version", stdout=TRUE, stderr=FALSE)[1]
  versionStr <- sub(".*version", "", versionStr)
  versionStr <- sub("^[^0-9]+", "", versionStr)
  versionPv <- try(package_version(versionStr))
  if (inherits(versionPv, "try-error")) {
    throw(sprintf("Unable to parse the '%s' version string, sorry!", command))
  } else {
    if (versionPv < version ||
        versionPv >= (version+1) )
      throw(sprintf("Failed to locate '%s', version '%s'.", command, version))
  }

  verbose && exit(verbose);

  pathname;
}) # findCmdv()


############################################################################
# HISTORY:
# 2013-03-07
# o Created TAT
############################################################################
