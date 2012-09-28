###########################################################################/**
# @RdocDefault systemGATK
#
# @title "Calls the GATK executable"
#
# \description{
#  @get "title".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments specifying GATK command line switches.}
#   \item{.fake}{If @TRUE, the executable is not called.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{\dontrun{
#   pathnameBAM <- "bwaData/LambdaVirusExample,bwa,is/Generic/reads_1.bam"
#   res <- systemGATK(T="CountReads", ..., stderr=FALSE)
# }}
#
#
# @author
#*/###########################################################################
setMethodS3("systemGATK", "default", function(..., .fake=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments '...':
  args <- list(...);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Calling GATK executable");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathname <- findGATK(verbose=less(verbose, 50));
  verbose && cat(verbose, "GATK jar file: ", pathname);
  pathname <- Arguments$getReadablePathname(pathname);

  # The actual binary is 'java'
  bin <- "java";


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Split up '...' arguments by system2() and GATK executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  keep <- is.element(names(args), names(formals(base::system2)));
  system2Args <- args[keep];
  args <- args[!keep];

  verbose && cat(verbose, "Arguments passed to system2():");
  verbose && str(verbose, system2Args);

  verbose && cat(verbose, "Arguments passed to java (to launch GATK):");
  
  args <- c(list(jar=sprintf("\"%s\"", pathname)), args);
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
}) # systemGATK()


############################################################################
# HISTORY:
# 2012-09-28
# o Created from systemPicard.R.
############################################################################
