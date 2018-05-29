###########################################################################/**
# @RdocDefault systemHTSeqCount
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
# \references{
#  [1] HTSeq: Analysing high-throughput sequencing data with Python,
#      June 2013.
#      \url{http://www-huber.embl.de/users/anders/HTSeq/}
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemHTSeqCount", "default", function(..., args=NULL, stdout=TRUE, stderr=FALSE, command="htseq-count", .fake=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'args':
  if (is.list(args)) {
    args <- unlist(args, use.names=TRUE);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling htseq-count executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locates the htseq-count executable
  bin <- findHTSeq(command=command, verbose=less(verbose, 50));
  verbose && cat(verbose, "Executable:");
  verbose && print(verbose, bin);

  verbose && cat(verbose, "Arguments passed to system2():");
  verbose && str(verbose, list(stdout=stdout, stderr=stderr));
  args <- c(list(...), args);
  verbose && cat(verbose, "Arguments passed to htseq-count:");
  verbose && print(verbose, args);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # System call
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Arguments as command-line options:");

  verbose && enter(verbose, "Calling system2()");
  cmdArgs <- sprintf("%s %s", names(args), args);
  cmdArgs <- trim(cmdArgs);
  cmdArgs <- cmdArgs[nchar(cmdArgs) > 0L];
  verbose && cat(verbose, "Arguments:");
  verbose && print(verbose, cmdArgs);
  if (!.fake) {
    res <- system2(bin, args=cmdArgs, stdout=stdout, stderr=stderr);
  } else {
    cat("<fake run>\n")
    res <- "<fake run>";
  }
  verbose && exit(verbose);

  verbose && exit(verbose);
  res;
}) # systemHTSeqCount()
