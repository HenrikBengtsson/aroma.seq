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
#   \item{...}{Additional arguments specifying Cuffdiff command line switches.}
#   \item{.fake}{If @TRUE, the executable is not called.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author
#*/###########################################################################

## cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]

setMethodS3("systemCuffdiff", "default", function(commandName="cuffdiff",
                                                  gtfFile,
                                                  samFiles,  ## vector with at least two sam files
                                                  ...,
                                                  system2ArgsList=list(stdout=FALSE),  ## For now, explicitly split off arguments to be passed to system2
                                                  .fake=FALSE, verbose=FALSE) {

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Arguments '...':
  dotArgs <- list(...);  ## i.e. dotArgs = list with one item named 'args'; these are the arguments to *cuffdiff*

  # Argument 'samFiles':
  if (length(samFiles) < 2) {
    throw("samFiles argument too short; need at least two sam files")
  }

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
  verbose && cat(verbose, "Executable: ", bin);

  verbose && cat(verbose, "Arguments passed to system2():");
  verbose && str(verbose, system2ArgsList)
  verbose && cat(verbose, "Arguments passed to Cuffdiff:");
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

## Usage:   cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]

############################################################################
# HISTORY:
# 2013-02-15
# o Created. TAT
############################################################################
