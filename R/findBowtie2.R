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

  versionPattern <- c("-version"=".*version ([0-9.]+).*");
  pathname <- findExternal(mustExist=mustExists, command=command, ..., versionPattern=versionPattern, verbose=verbose);

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
