###########################################################################/**
# @RdocDefault systemBowtie2
#
# @title "Wrapper for bowtie2"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{commandName}{A @character string specifying the bowtie2 command.}
#   \item{...}{Additional arguments specifying bowtie2 command line switches.}
#   \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
#   \item{.fake}{If @TRUE, the executable is not called.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemBowtie2", "default", function(commandName="bowtie2",
                                                 ...,
                                                 system2ArgsList=list(stdout=TRUE, stderr=FALSE),
                                                 .fake=FALSE, verbose=FALSE
                                                 ) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments '...':
  dotArgs <- list(...);  ## list with one item named 'args'; these are the arguments to *bowtie2*

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling bowtie2 executable");

  ## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate the bowtie2 executable
  bin <- findBowtie2(command=commandName, verbose=less(verbose, 50));
  verbose && cat(verbose, "Executable: ", bin);
  verbose && cat(verbose, "Arguments passed to system2():");
  verbose && str(verbose, system2ArgsList)
  verbose && cat(verbose, "Arguments passed to bowtie2:");
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
  callArgs <- c(callArgs, system2ArgsList);

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

  ## Strategy:  Previously, built up the command like this
  ##  cmdOptions <- .optionsList2String(optionsList)
  ##  cmdOptionsAdd <- .optionsList2String(list(...))
  ##  cmdOptions <- paste(cmdOptions, cmdOptionsAdd)
  ## Now assume the argument string is passed fully formed to systemBowtie2.

})
