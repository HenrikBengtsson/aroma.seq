## source("systemTophatBuild.R")
## - 20130215 This code will likely be superseded by systemTophat.R; call that from a 'tophatBuild.R' (to be written)

if (FALSE) {  ## (for standalone testing)
  library(R.filesets)
  source("findTopHat.R")
}

###########################################################################/**
# @RdocDefault systemTopHatBuild
#
# @title "Calls the TopHat executable, specifically to build a transcriptome index"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments specifying TopHat command line switches.}
#   \item{.fake}{If @TRUE, the executable is not called.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author
#*/###########################################################################
setMethodS3("systemTopHatBuild", "default", function( ## ( No 'command' arg in this case )
                                                     ...,
                                                     system2Args,  ## explicitly split off arguments to be passed to system2
                                                     .fake=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Arguments '...':
  args <- list(...);  ## These are the arguments to *tophat*

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling TopHat executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bin <- findTopHat(verbose=less(verbose, 50));
  verbose && cat(verbose, "Executable: ", bin);

  verbose && cat(verbose, "Arguments passed to system2():");
  verbose && str(verbose, system2Args);
  verbose && cat(verbose, "Arguments passed to TopHat:");
  verbose && str(verbose, args);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup command line switches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## args <- paste(keys, args, sep=" ");
  args <- trim(args);
  verbose && cat(verbose, "Command line options:");
  verbose && print(verbose, args);

  nms <- optionsList
  nms <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")
  names(optionsList) <- nms

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # System call
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "System call:");
  cmd <- sprintf("%s %s", bin, paste(args, collapse=" "));
  verbose && print(verbose, cmd);
  verbose && str(verbose, system2Args);

  verbose && enter(verbose, "system2() call");
  callArgs <- list(command=bin, args=args);
  callArgs <- c(callArgs, system2Args);
  verbose && str(verbose, callArgs);
  if (!.fake) {
    res <- do.call(base::system2, callArgs);
  } else {
    res <- "<fake run>";
  }

  verbose && exit(verbose);
  verbose && exit(verbose);
  res;
}) # systemTopHatBuild()


############################################################################
# HISTORY:
# 2013-01-24
# o Created. TAT
############################################################################
