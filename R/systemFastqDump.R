###########################################################################/**
# @RdocDefault systemFastqDump
#
# @title "Wrapper for fastq-dump SRA utility"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{commandName}{A @character string specifying the FastqDump command.}
#   \item{...}{Additional arguments specifying FastqDump command line switches.}
#   \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
#   \item{.fake}{If @TRUE, the executable is not called.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemFastqDump", "default", function(...,
                                                   commandName="fastq-dump",
                                                   Stdout=TRUE,
                                                   Stderr=FALSE,
                                                   .fake=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # Arguments '...':
  Args <- unlist(list(...))
  # Args <- Args$args; # char vector of args passed to fastq-dump

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling FastqDump executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bin <- findFastqDump(commandName=commandName);
  verbose && cat(verbose, "Executable:");
  verbose && print(verbose, bin);

  # verbose && cat(verbose, "Arguments passed to system2():");
  # verbose && str(verbose, system2ArgsList)
  verbose && cat(verbose, "Arguments passed to FastqDump:");
  verbose && str(verbose, Args);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # System call
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  verbose && enter(verbose, "system2() call");
  cmd <- sprintf("%s %s", bin, paste(Args, collapse=" "));
  verbose && print(verbose, cmd);

  if (!.fake) {
    res <- system2(command=commandName, args=Args, stderr=Stderr, stdout=Stdout)
  } else {
    cat("<fake run>\n")
    res <- "<fake run>";
  }

  verbose && exit(verbose);
  verbose && exit(verbose);
  res;
}) # systemFastqDump()
