###########################################################################/**
# @set class=GenericDataFileSet
# @RdocMethod dsApply
#
# @title "Applies a function to each file in the file set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{FUN}{A @function.}
#  \item{...}{Arguments passed to \code{FUN}.}
#  \item{args}{(optional) A @list of additional arguments
#    passed to \code{FUN}.}
#  \item{skip}{If @TRUE, already processed files are skipped.}
#  \item{verbose}{See @see "R.utils::Verbose".}
#  \item{.parallel}{If @TRUE, parallel/distributed processing is utilized.}
#  \item{.control}{(internal) A named @list structure controlling
#        the processing.}
# }
#
# \value{
#   Returns nothing.
# }
#
# \seealso{
#   The \pkg{BatchJobs} package is utilized for parallel/distributed
#   processing.
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("dsApply", "GenericDataFileSet", function(ds, FUN, ..., args=list(), skip=FALSE, verbose=FALSE, .parallel=FALSE, .control=list(dW=1.0)) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'FUN':
  stopifnot(is.function(FUN));

  # Arguments '...':
  vargs <- list(...);
  nvargs <- length(vargs);

  # Argument 'args':
  if (!is.list(args)) {
    throw("Argument 'args' must be a list: ", mode(args));
  }

  # Argument 'skip':
  skip <- Arguments$getLogical(skip);

  # Argument '.parallel':
  .parallel <- Arguments$getLogical(.parallel);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }



  # Run via BatchJobs?
  parallel <- getOption(aromaSettings, "devel/BatchJobs", .parallel);

  # The additional set of arguments passed in each function call
  vargs <- c(vargs, args);
  allArgs <- c(vargs, skip=skip, verbose=verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Alt 1: BatchJobs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (parallel) {
    verbose && enter(verbose, "Processing using BatchJobs");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Setup
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Poll roughly every dW seconds
    dW <- .control$dW;
    if (is.null(dW)) dW <- 1.00;
    dW <- Arguments$getNumeric(dW, range=c(0,Inf));

    # BatchJob registry to be used
    reg <- .getBatchJobRegistry(ds, args=vargs);
    verbose && print(verbose, reg);


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (i) Add jobs, iff missing
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nbrOfJobs <- getJobNr(reg);
    if (nbrOfJobs == 0L) {
      verbose && enter(verbose, "Adding jobs to registry");
      # Tweak FUN()
##      FUNx <- function(...) {
##        # Allow comma-separated registry directory paths
##        oopts <- options("BatchJobs.check.posix");
##        on.exit({ options(oopts) }, add=TRUE);
##        options("BatchJobs.check.posix"=FALSE);
##        FUN(...);
##      } # FUNx()
      ids <- batchMap(reg, fun=FUN, getFiles(ds), more.args=allArgs);
      verbose && cat(verbose, "Job IDs added:");
      verbose && str(verbose, ids);
      verbose && print(verbose, reg);
      verbose && exit(verbose);
    }

    verbose && print(verbose, showStatus(reg));
#    throw("Jobs have already been added: ", reg$id);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (ii) Launch jobs
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Launching jobs");
    lastTodo <- NULL;
    todo <- findNotSubmitted(reg);
    if (length(todo) > 0L) {
      # (a) Wait and see if jobs are being submitted by other process
      while(!identical(todo, lastTodo)) {
         lastTodo <- todo;
         Sys.sleep(1.0);
         todo <- findNotRunning(reg);
      }

      verbose && cat(verbose, "Job IDs to be submitted:");
      verbose && print(verbose, todo);
      submitted <- submitJobs(reg, ids=todo);
      verbose && cat(verbose, "Job IDs actually submitted:");
      verbose && print(verbose, submitted);
      verbose && cat(verbose, "Job IDs not submitted:");
      verbose && print(verbose, setdiff(todo, submitted));
    } else {
      verbose && cat(verbose, "No new jobs to be submitted.");
    }
    verbose && exit(verbose);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (iii) Wait for jobs to finish
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Waiting for jobs to finish");
    t0 <- Sys.time();
    tCount <- 0L;
    status <- NULL;
    while (length(findNotTerminated(reg)) > 0L) {
      lastStatus <- status;
      out <- capture.output(status <- showStatus(reg));
      if (identical(status, lastStatus)) {
        verbose && writeRaw(verbose, ".");
        # Time stamp?
        dt <- difftime(Sys.time(), t0, units="secs");
        dMins <- as.integer(dt) %/% 10;
        if (dMins > tCount) {
          tCount <- dMins;
          if (dt > 1.5*60) {
            units(dt) <- "mins";
          } else if (dt > 1.5*3600) {
            units(dt) <- "hours";
          }
          verbose && writeRaw(verbose, sprintf("[%s]\n", format(dt)));
        }
      } else {
        verbose && writeRaw(verbose, "\n");
        verbose && print(verbose, status);
      }
      Sys.sleep(dW);
    } # while(...)
    verbose && exit(verbose);

    verbose && print(verbose, showStatus(reg));

    verbose && exit(verbose);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Alt 2: Sequentially using regular R
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!parallel) {
    for (ii in seq_along(ds)) {
      df <- getFile(ds, ii);
      verbose && enter(verbose, sprintf("Item #%d ('%s') of %d", ii, getName(df), length(ds)));
      argsII <- c(list(df), allArgs);
      verbose && cat(verbose, "Call arguments:");
      verbose && str(verbose, argsII);
      res <- do.call(FUN, args=argsII);
      verbose && str(verbose, res);

      # Not needed anymore
      df <- argsII <- res <- NULL;

      verbose && exit(verbose);
    } # for (ii ...)
  }

  res <- NULL;
  invisible(res);
}, protected=TRUE) # dsApply()


setMethodS3(".getBatchJobRegistryId", "default", function(class, label=NULL, version=NULL, ..., verbose=FALSE) {
  # Argument 'class':
  class <- Arguments$getCharacters(class);

  # Argument 'label':
  label <- Arguments$getCharacter(label);

  # Argument 'version':
  version <- Arguments$getCharacter(version);


  # Construct key from all object/arguments
  key <- list(
    class=class,
    label=label,
    version=version,
    ...
  );
  keyId <- R.cache::getChecksum(key, algo="crc32");
  keyTime <- format(Sys.time(), "%Y%m%d%H%M%S");
  pid <- Sys.getpid();
  id <- paste(c(key$class[1L], key$label, keyId, keyTime, pid), collapse="_");
  id <- Arguments$getCharacter(id, length=c(1L,1L));

  id;
}, protected=TRUE)


setMethodS3(".getBatchJobRegistryId", "GenericDataFileSet", function(object, ...) {
  files <- getFiles(object);
  keys <- lapply(files, FUN=function(file) {
    list(filename=getFilename(file), fileSize=getFileSize(file));
  });
  .getBatchJobRegistryId(class=class(object), fileKeys=keys, ...);
}, protected=TRUE)



setMethodS3(".getBatchJobRegistry", "default", function(..., skip=TRUE) {
  .useBatchJobs();
  require("BatchJobs"); # To please R CMD check.

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'skip':
  skip <- Arguments$getLogical(skip);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Constructor BatchJobs registry ID...
  id <- .getBatchJobRegistryId(...);

  # ...and registry directory
  rootPath <- ".batchJobsRegistries";
  path <- file.path(rootPath, id);
  path <- Arguments$getWritablePath(path);

  # Allow comma-separated registry directory paths
  oopts <- options("BatchJobs.check.posix");
  on.exit({ options(oopts) }, add=TRUE);
  options("BatchJobs.check.posix"=FALSE);

  # Load BatchJobs registry or create if missing
  packages <- c("aroma.seq");
  reg <- BatchJobs::makeRegistry(id=id, file.dir=path, skip=skip, packages=packages);

  reg;
}, protected=TRUE) # .getBatchJobRegistry()


.useBatchJobs <- function(version="*", ...) {
  if (identical(version, "*")) {
    desc <- packageDescription("aroma.seq");
    desc <- desc[c("Depends", "Imports", "Suggests")];
    desc <- unlist(desc, use.names=FALSE);
    desc <- strsplit(desc, split=",", fixed=TRUE);
    desc <- unlist(desc, use.names=FALSE);
    desc <- gsub("\n", "", desc, fixed=TRUE);
    pattern <- "[ ]*([^ (]+)[ ]*(()|[(][^]]+[)])";
    pkgs <- gsub(pattern, "\\1", desc);
    vers <- gsub(pattern, "\\2", desc);
    vers <- gsub("[()]", "", vers);
    names(vers) <- pkgs;
    dups <- duplicated(pkgs);
    vers <- vers[!dups];
    version <- vers["BatchJobs"];
  }
  R.utils::use("BatchJobs", version=version, ...);
} # .useBatchJobs()



############################################################################
# HISTORY:
# 2013-09-28
# o ROBUSTNESS: Now .getBatchJobRegistryId() adds process ID ("pid") and#   a timestamp to the registry path to make it even more unique.
# o BUG FIX: dsApply() for GenericDataFileSet when executed via the
#   'BatchJobs' methods would not allow commas in the work directory
#   among other directories.
# 2013-08-31
# o Now .useBatchJobs() utilizes R.utils::use().
# o Added dsApply() for GenericDataFileSet.  Added Rdocs.
# 2013-08-26
# o Added .getBatchJobRegistry() and .getBatchJobRegistryId().
# o Added .usePackage() and .useBatchJobs().
# o Created.
############################################################################
