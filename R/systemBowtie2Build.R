###########################################################################/**
# @RdocDefault systemBowtie2Build
#
# @title "Wrapper for bowtie2-build"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{commandName}{A @character string specifying the bowtie2 build command.}
#   \item{...}{Additional arguments specifying bowtie2-build command line switches.}
#   \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
#   \item{.fake}{If @TRUE, the executable is not called.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemBowtie2Build", "default", function(commandName="bowtie2-build",
                                                      ...,
                                                      system2ArgsList=list(stdout=TRUE, stderr=FALSE),
                                                      .fake=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments '...':
  dotArgs <- list(...);  ## list with one item named 'args'; these are the arguments to *bowtie2-build*

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling bowtie2-build executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locates the bowtie2-build executable
  bin <- findBowtie2(command=commandName, verbose=less(verbose, 50));
  verbose && cat(verbose, "Executable:");
  verbose && print(verbose, bin);
  verbose && cat(verbose, "Arguments passed to system2():");
  verbose && str(verbose, system2ArgsList)
  verbose && cat(verbose, "Arguments passed to bowtie2-build:");
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
}) # systemBowtie2Build()
