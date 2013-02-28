###########################################################################/**
# @RdocDefault systemTopHat
#
# @title "Calls the TopHat executable to align input reads"
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
setMethodS3("systemTopHat", "default", function(commandName,
                                                ...,
                                                system2ArgsList=list(),  ## For now, explicitly split off arguments to be passed to system2
                                                .fake=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Arguments '...':
  dotArgs <- list(...);  ## list with one item named 'args'; these are the arguments to *tophat*

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
  bin <- findCmd(commandName, verbose=less(verbose, 50));
  verbose && cat(verbose, "Executable: ", bin);

  verbose && cat(verbose, "Arguments passed to system2():");
  verbose && str(verbose, system2ArgsList)
  verbose && cat(verbose, "Arguments passed to TopHat:");
  verbose && str(verbose, dotArgs);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup command line switches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  ## dotArgs <- trim(dotArgs);
  verbose && cat(verbose, "Command line options:");
  verbose && print(verbose, dotArgs);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # System call
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "System call:");
  cmd <- sprintf("%s %s", bin, paste(dotArgs, collapse=" "));
  verbose && print(verbose, cmd);
  verbose && str(verbose, system2ArgsList);

  verbose && enter(verbose, "system2() call");
  callArgs <- list(command=bin, args=paste(names(dotArgs$args), dotArgs$args, sep=" "))

  verbose && str(verbose, callArgs);
  if (!.fake) {
    res <- do.call(what=base::system2, args=callArgs);
  } else {
    res <- "<fake run>";
  }

  verbose && exit(verbose);
  verbose && exit(verbose);
  res;
}) # systemTopHat()


############################################################################
# HISTORY:
# 2013-01-24
# o Created. TAT
############################################################################
