###########################################################################/**
# @RdocDefault systemCuffdiff
#
# @title "Calls the Cuffdiff executable to perform isoform abundance estimation"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{cuffdiffStr}{A @character string of arguments to the cuffdiff command}
#   \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
#   \item{...}{...}
#   \item{commandName}{Name of executable}
#   \item{.fake}{If @TRUE, the executable is not called.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################

## cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]

setMethodS3("systemCuffdiff", "default", function(cuffdiffStr="",
                                                  system2ArgsList=list(stdout=FALSE),
                                                  ...,
                                                  commandName="cuffdiff",
                                                  .fake=TRUE, verbose=TRUE) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Arguments 'commandName':
  if (is.null(commandName)) {
    throw("Argument 'commandName' is null; supply the name of a command to run")
  }

  # Arguments '...':
  dotArgs <- list(...);  ## [ Not used ]

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling Cuffdiff executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bin <- findCmd(commandName, verbose=less(verbose, 50));
  verbose && cat(verbose, "Executable:");
  verbose && print(verbose, bin);
  verbose && cat(verbose, "Arguments passed to system2():");
  verbose && str(verbose, system2ArgsList)
  verbose && cat(verbose, "Arguments passed to Cuffdiff:");
  verbose && str(verbose, cuffdiffStr)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # System call
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Command line args:");
  verbose && print(verbose, cuffdiffStr);
  verbose && cat(verbose, "System call:");
  cmd <- sprintf("%s %s", bin, cuffdiffStr)  ## << this is the command line
  verbose && print(verbose, cmd);

  verbose && enter(verbose, "system2() call");
  callArgs <- list(command=bin, args=cuffdiffStr)
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
}) # systemCuffdiff()
