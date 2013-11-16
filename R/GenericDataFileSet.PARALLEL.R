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
#  \item{.parallel}{A @character string specifying what mechanism to use
#    for performing parallel processing, if at all.}
#  \item{.control}{(internal) A named @list structure controlling
#        the processing.}
# }
#
# \value{
#   Returns nothing.
# }
#
# \seealso{
#   The \pkg{BiocParallel} and \pkg{BatchJobs} packages is utilized
#   for parallel/distributed processing, depending on settings.
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("dsApply", "GenericDataFileSet", function(ds, FUN, ..., args=list(), skip=FALSE, verbose=FALSE, .parallel=c("none", "BatchJobs", "BiocParallel::BatchJobs"), .control=list(dW=1.0)) {
  # To please R CMD check (because BatchJobs is just "suggested")
  getJobNr <- batchMap <- showStatus <- findNotSubmitted <-
      findNotRunning <- submitJobs <- findNotTerminated <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  assertNoGlobalVariables <- function(FUN, ...) {
    # TO DO...
    ## globals <- findGlobals(FUN, merge=FALSE);
  } # assertNoGlobalVariables()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'FUN':
  stopifnot(is.function(FUN));
  assertNoGlobalVariables(FUN);


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
  parallel <- getOption(aromaSettings, "devel/parallel", .parallel);
  parallel <- match.arg(parallel, choices=eval(formals(dsApply.GenericDataFileSet)$.parallel));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Processing ", class(ds)[1L]);
  verbose && cat(verbose, "Mechanism for parallel processing: ", parallel);
  verbose && print(verbose, ds);

  # The additional set of arguments passed in each function call
  vargs <- c(vargs, args);
  allArgs <- c(vargs, skip=skip, verbose=verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Alt 1: Sequentially using regular R
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- NULL;
  if (parallel == "none") {
    for (ii in seq_along(ds)) {
      df <- ds[[ii]];
      verbose && enter(verbose, sprintf("Item #%d ('%s') of %d", ii, getName(df), length(ds)));

      argsII <- c(list(df), allArgs);
      verbose && cat(verbose, "Call arguments:");
      verbose && str(verbose, argsII);
      argsII$verbose <- less(verbose, 1);
      res <- do.call(FUN, args=argsII);
      verbose && str(verbose, res);

      # Not needed anymore
      df <- argsII <- res <- NULL;

      verbose && exit(verbose);
    } # for (ii ...)
  } # if (parallel == "none")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Alt 2: BatchJobs processing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (parallel == "BatchJobs") {
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
      ids <- batchMap(reg, fun=FUN, as.list(ds), more.args=allArgs);
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
  } # if (parallel == "BatchJobs")


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Alt 3: BiocParallel processing
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (parallel == "BiocParallel::BatchJobs") {
    verbose && enter(verbose, "Processing using BiocParallel");

    # WORKAROUND: Make sure 'methods' package is *attached*, not
    # just loaded. /HB 2013-11-09
    pkgName <- "methods";
    require(pkgName, character.only=TRUE) || throw("Package not attached: ", pkgName);

    # WORKAROUND: Make sure to load BatchJobs config files
    require("BatchJobs", character.only=TRUE) || throw("Package not attached: BatchJobs");

    conffile <- c(".BatchJobs.R", "~/.BatchJobs.R")
    conffile <- normalizePath(conffile, mustWork=FALSE);
    conffile <- conffile[file_test("-f", conffile)];
    if (length(conffile) > 0L) {
      conffile <- conffile[1L];
      if (isFile(conffile)) loadConfig(conffile);
    }

    # WORKAROUND: Allow for commas in BatchJobs-related pathnames
    oopts <- options("BatchJobs.check.posix");
    on.exit({ options(oopts) }, add=TRUE);
    options("BatchJobs.check.posix"=FALSE);

    bpParam <- BatchJobsParam();
    register(bpParam);
    verbose && cat(verbose, "Using parameters:");
    verbose && print(verbose, bpParam);

    verbose && enter(verbose, "Calling bplapply()");
    args <- c(list(ds, FUN=FUN), allArgs, BPPARAM=bpParam);
    verbose && cat(verbose, "Arguments passed to bplapply():");
    verbose && str(verbose, args);
    res <- do.call(BiocParallel::bplapply, args);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # if (parallel == "BiocParallel")

  verbose && exit(verbose);

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
  keys <- lapply(object, FUN=function(file) {
    list(filename=getFilename(file), fileSize=getFileSize(file));
  });
  .getBatchJobRegistryId(class=class(object), fileKeys=keys, ...);
}, protected=TRUE)



setMethodS3(".getBatchJobRegistry", "default", function(..., skip=TRUE) {
  .useBatchJobs();
  # To please R CMD check.
  require("BatchJobs");

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
# 2013-11-02
# o Adding support for distributed processing via 'BiocParallel'.
# 2013-11-01
# o Made the code more generic.
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
