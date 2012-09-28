###########################################################################/**
# @RdocDefault systemPicard
#
# @title "Calls the Picard executable"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#   \item{command}{A @character string specifying the Picard command
#     (the name of the jar file without the *.jar extension).}
#   \item{...}{Additional arguments specifying Picard command line switches.}
#   \item{.fake}{If @TRUE, the executable is not called.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author
#*/###########################################################################
setMethodS3("systemPicard", "default", function(command, ..., .fake=FALSE, verbose=FALSE) {
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

  verbose && enter(verbose, "Calling Picard executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- findPicard(verbose=less(verbose, 50));
  verbose && cat(verbose, "Picard directory: ", path);

  verbose && cat(verbose, "Picard command: ", command);
  filename <- sprintf("%s.jar", command);

  pathname <- file.path(path, filename);
  verbose && cat(verbose, "Pathname to Picard jar file: ", pathname);
  pathname <- Arguments$getReadablePathname(pathname);

  # The actual binary is 'java'
  bin <- "java";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Split up '...' arguments by system2() and Picard executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  keep <- is.element(names(args), names(formals(base::system2)));
  system2Args <- args[keep];
  args <- args[!keep];

  verbose && cat(verbose, "Arguments passed to system2():");
  verbose && str(verbose, system2Args);

  verbose && cat(verbose, "Arguments passed to java (to launch Picard):");
  
  args <- c(list(sprintf("-jar \"%s\"", pathname)), args);
  verbose && str(verbose, args);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup command line switches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Automatically add dashes
  keys <- names(args);
  bad <- grep("^-", keys, value=TRUE);
  if (length(bad) > 0L) {
    throw("Detected non-valid command line options: ", hpaste(bad));
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
}) # systemPicard()


############################################################################
# HISTORY:
# 2012-09-27
# o Created.
############################################################################
