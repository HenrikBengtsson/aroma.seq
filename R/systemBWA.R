###########################################################################/**
# @RdocDefault systemBWA
#
# @title "Calls the BWA executable"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{command}{A @character string specifying the BWA command.}
#   \item{...}{Additional arguments specifying BWA command line switches.}
#   \item{.fake}{If @TRUE, the executable is not called.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemBWA", "default", function(command, ..., .fake=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'command':
  command <- Arguments$getCharacter(command);

  # Arguments '...':
  args <- list(...);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling BWA executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bin <- findBWA(verbose=less(verbose, 50));
  verbose && cat(verbose, "Executable: ", bin);

  # Unwrap vectors
  args <- unlist(args);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Split up '...' arguments by system2() and BWA executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  keys <- names(args);
  if (is.null(keys)) {
    system2Args <- NULL;
  } else {
    keep <- is.element(keys, names(formals(base::system2)));
    system2Args <- args[keep];
    args <- args[!keep];
  }

  verbose && cat(verbose, "Arguments passed to system2():");
  verbose && str(verbose, system2Args);

  verbose && cat(verbose, "Arguments passed to BWA:");
  args <- c(list(command), args);
  verbose && str(verbose, args);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup command line switches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Automatically add dashes
  keys <- names(args);
  missing <- grep("^[^-]", keys);
  if (length(missing) > 0L) {
    invalid <- (args[missing] == "");
    if (any(invalid)) {
      throw("Detected non-valid command line switched: ", hpaste(args[missing][invalid]));
    }
    keys[missing] <- sprintf("-%s", keys[missing]);
  }

  args <- paste(keys, args, sep=" ");
  args <- trim(args);
  verbose && cat(verbose, "Command line options:");
  verbose && print(verbose, args);


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
}) # systemBWA()
