## Create a 'systemHtseqCount.R", modeled on systemBowtie2Build.R


###########################################################################/**
# @RdocDefault systemHtseqCount
#
# @title "Wrapper for htseq-count"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{commandName}{A @character string specifying the htseq-count command.}
#   \item{...}{Additional arguments specifying htseq-count command line switches.}
#   \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
#   \item{.fake}{If @TRUE, the executable is not called.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemHtseqCount", "default", function(commandName="htseq-count",
                                                    ...,
                                                    system2ArgsList=list(stdout=TRUE, stderr=FALSE),
                                                    .fake=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments '...':
  dotArgs <- list(...);  ## list with one item named 'args'; these are the arguments to *htseq-count*

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling htseq-count executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locates the htseq-count executable
  bin <- findHtseq(command=commandName, verbose=less(verbose, 50));
  verbose && cat(verbose, "Executable: ", bin);
  verbose && cat(verbose, "Arguments passed to system2():");
  verbose && str(verbose, system2ArgsList)
  verbose && cat(verbose, "Arguments passed to htseq-count:");
  verbose && str(verbose, dotArgs);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup command line switches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
  callArgs <- c(callArgs, system2ArgsList)

  verbose && str(verbose, callArgs);
  if (!.fake) {
    res <- do.call(what=base::system2, args=callArgs);
  } else {
    cat("<fake run>\n")
    res <- "<fake run>";
  }

  verbose && exit(verbose);
  verbose && exit(verbose);
  res;
}) # systemHtseqCount()


############################################################################
# HISTORY:
# 2013-05-31
# o TAT:  Created systemHtseqCount stub
############################################################################
