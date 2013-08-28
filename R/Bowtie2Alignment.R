###########################################################################/**
# @RdocClass Bowtie2Alignment
#
# @title "The Bowtie2Alignment class"
#
# \description{
#  @classhierarchy
#
#  ...
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "AbstractAlignment".}
#  \item{indexSet}{An @see "Bowtie2IndexSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Supported operating systems}{
#   This method is available on Linux, OSX, and Windows [1].
# }
#
# \author{Henrik Bengtsson and Pierre Neuvial}
#
# \references{
#  [1] Bowtie2, John Hopkins University, 2013.
#      \url{http://bowtie-bio.sourceforge.net/bowtie2/}
# }
#*/###########################################################################
setConstructorS3("Bowtie2Alignment", function(..., indexSet=NULL) {
  # Validate arguments
  if (!is.null(indexSet)) {
    indexSet <- Arguments$getInstanceOf(indexSet, "Bowtie2IndexSet");
  }

  # Arguments '...':
  args <- list(...);

  extend(AbstractAlignment(..., indexSet=indexSet), "Bowtie2Alignment");
})


setMethodS3("getParameters", "Bowtie2Alignment", function(this, ...) {
  params <- NextMethod("getParameters");
  params <- c(params, getOptionalArguments(this, ...));
  params;
}, protected=TRUE)


###########################################################################/**
# @RdocMethod process
#
# @title "Runs the aligner"
#
# \description{
#   @get "title" on all input files.
#   The generated BAM files are all sorted and indexed.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
#  \item{skip}{If @TRUE, already processed files are skipped.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "BamDataSet".
# }
#
# @author
#*/###########################################################################
setMethodS3("process", "Bowtie2Alignment", function(this, ..., skip=TRUE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  asBowtie2Parameters <- function(rg, ...) {
    if (isEmpty(rg)) {
      return(NULL);
    }

    # Validate
##    if (!hasID(rg)) {
##      throw("Bowtie requires that the SAM read group has an ID.");
##    }

    rgId <- asSamList(rg)$ID;
    if (is.null(rgId)) rgId <- 1L;

    rgArgs <- asString(rg, fmtstr="%s:%s");
    rgArgs <- rgArgs[regexpr("^ID:", rgArgs) == -1L];

    # Don't forget to put within quotation marks
    rgArgs <- sprintf("\"%s\"", rgArgs);

    rgArgs <- as.list(rgArgs);
    names(rgArgs) <- rep("--rg", times=length(rgArgs));

    rgArgs <- c(list("--rg-id"=rgId), rgArgs);

    rgArgs;
  } # asBowtie2Parameters()


  bowtie2 <- function(pathnameFQ, indexPrefix, pathnameSAM, ..., verbose=FALSE) {
    pathnameFQ <- Arguments$getReadablePathnames(pathnameFQ, length=1:2);
    isPaired <- (length(pathnameFQ) == 2L);

    indexPrefix <- Arguments$getCharacter(indexPrefix);
    indexPath <- Arguments$getReadablePath(getParent(indexPrefix));
    pathnameSAM <- Arguments$getWritablePathname(pathnameSAM);

    # Check if FASTQ.gz files are supported
    isGzipped <- any(sapply(pathnameFQ, FUN=isGzipped));
    if (isGzipped) {
      if (is.na(gzAllowed)) {
        gzAllowed <- queryBowtie2("support:fastq.gz");
      }

      if (!gzAllowed) {
        decompress <- getOption(aromaSettings, "devel/fastq.gz/decompress", TRUE);
        if (!decompress) {
          why <- attr(gzAllowed, "why");
          throw(sprintf("Cannot align reads in '%s': %s", getPathname(df), why));
        }

        # If not, temporarily decompress (=remove when done)
        pathnameFQ <- sapply(pathnameFQ, FUN=gunzip, temporary=TRUE, remove=FALSE);
        on.exit({
          # Make sure to remove temporary file
          lapply(pathnameFQ, FUN=function(pathname) {
            if (isFile(pathname)) file.remove(pathname);
          });
        }, add=TRUE);

        # Sanity check
        isGzipped <- any(sapply(pathnameFQ, FUN=isGzipped));
        stopifnot(!isGzipped);
      }
    } # if (isGzipped)

    # WORKAROUND: Bowtie2() does not support commas in the FASTQ
    # pathname.  If so, use a temporary filename without commas.
    hasComma <- (regexpr(",", pathnameFQ, fixed=TRUE) != -1L);
    if (any(hasComma)) {
      pathnameFQ[hasComma] <- sapply(pathnameFQ[hasComma], FUN=function(pathname) {
        ext <- if (isGzipped(pathname)) ".fq.gz" else ".fq";
        pathnameT <- tempfile(fileext=ext);
        createLink(pathnameT, target=pathname);
        pathnameT;
      });

      # Remove temporary files
      on.exit({
        file.remove(pathnameFQ[hasComma]);
      }, add=TRUE);
    }

    if (isPaired) {
      # Sanity check
      stopifnot(pathnameFQ[1L] != pathnameFQ[2L]);
      res <- systemBowtie2(args=list("-x"=indexPrefix, "-1"=pathnameFQ[1L], "-2"=pathnameFQ[2L], "-S"=pathnameSAM, ...), verbose=verbose);
    } else {
      res <- systemBowtie2(args=list("-x"=indexPrefix, "-U"=pathnameFQ, "-S"=pathnameSAM, ...), verbose=verbose);
    }

    # In case bowtie2 generates empty SAM files
    # /HB 2012-10-01 (still to be observed)
    if (isFile(pathnameSAM)) {
      if (file.info(pathnameSAM)$size == 0L) {
        verbose && cat(verbose, "Removing empty SAM file falsely created by Bowtie2: ", pathnameSAM);
        file.remove(pathnameSAM);
      }
    }

    res;
  } # bowtie2()


  processOne <- function(df, paired=FALSE, indexPrefix, rgSet, params, ..., path, skip=TRUE) {
    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    verbose && enter(verbose, "Bowtie2 alignment of one sample");

    # Paired-end?
    if (paired) {
      pair <- list(R1=df, R2=getMateFile(df));
      pathnameFQ <- sapply(pair, FUN=getPathname);
      verbose && cat(verbose, "FASTQ R1 pathname: ", pathnameFQ[1L]);
      verbose && cat(verbose, "FASTQ R2 pathname: ", pathnameFQ[2L]);
    } else {
      pathnameFQ <- getPathname(df);
      verbose && cat(verbose, "FASTQ pathname: ", pathnameFQ);
    }

    # The SAM and BAM files to be generated
    fullname <- getFullName(df);
    filename <- sprintf("%s.sam", fullname);
    pathnameSAM <- Arguments$getWritablePathname(filename, path=path);
    verbose && cat(verbose, "SAM pathname: ", pathnameSAM);
    filename <- sprintf("%s.bam", fullname);
    pathnameBAM <- Arguments$getWritablePathname(filename, path=path);
    verbose && cat(verbose, "BAM pathname: ", pathnameBAM);

    done <- (skip && isFile(pathnameBAM));
    if (done) {
      verbose && cat(verbose, "Already aligned. Skipping");
    } else {
      verbose && print(verbose, df);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (a) Generate SAM file
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!isFile(pathnameSAM)) {
        # Extract sample-specific read group
        rgII <- getSamReadGroup(df);
        if (length(rgSet) > 0L) {
          rgII <- merge(rgSet, rgII);
        }
        verbose && cat(verbose, "Writing SAM Read Groups:");
        verbose && print(verbose, rgII);
        verbose && cat(verbose, "Bowtie2 parameter:");
        rgArgs <- asBowtie2Parameters(rgII);
        verbose && print(verbose, rgArgs);
        rgII <- NULL; # Not needed anymore

        args <- list(
          pathnameFQ,
          indexPrefix=indexPrefix,
          pathnameSAM=pathnameSAM
        );
        args <- c(args, rgArgs);
        args <- c(args, params);
        verbose && cat(verbose, "Arguments:");
        verbose && str(verbose, args);
        args$verbose <- less(verbose, 5);

        res <- do.call(bowtie2, args=args);
        verbose && cat(verbose, "System result code: ", res);
      } # if (!isFile(pathnameSAM))

      # Sanity check
      stopifnot(isFile(pathnameSAM));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (b) Generate BAM file from SAM file
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!isFile(pathnameBAM)) {
        sf <- SamDataFile(pathnameSAM);
        bf <- convertToBamDataFile(sf, verbose=less(verbose, 5));
        print(pathnameBAM);
      }
      # Sanity check
      stopifnot(isFile(pathnameBAM));
    }

    verbose && exit(verbose);

    invisible(list(pathnameFQ=pathnameFQ, pathnameSAM=pathnameSAM, pathnameBAM=pathnameBAM));
  } # processOne()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Bowtie2 alignment");

  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  is <- getIndexSet(this);
  verbose && cat(verbose, "Aligning using index set:");
  verbose && print(verbose, is);
  indexPrefix <- getIndexPrefix(is);

  rgSet <- this$.rgSet;
  if (!is.null(rgSet)) {
    verbose && cat(verbose, "Assigning SAM read group:");
    verbose && print(verbose, rgSet);
    validate(rgSet);
  }

  # Indicates whether gzipped FASTQ files are supported or not.
  gzAllowed <- NA;

  params <- getParameters(this);
  verbose && cat(verbose, "Additional bowtie2 arguments:");
  verbose && str(verbose, params);

  verbose && cat(verbose, "Number of files: ", length(ds));

  isPaired <- isPaired(this);
  if (isPaired) {
    pairs <- getFilePairs(ds);
  }

  outPath <- getPath(this);

  useBatchJobs <- getOption(aromaSettings, "devel/BatchJobs", FALSE);

  if (useBatchJobs) {
    verbose && enter(verbose, "Processing using BatchJobs");

    # BatchJob registry to be used
    dataSetArgs <- list(paired=isPaired, indexPrefix=indexPrefix, rgSet=rgSet, params=params, path=outPath);
    reg <- getBatchJobRegistry(ds, method="process", dataSetArgs=dataSetArgs);
    verbose && print(verbose, reg);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (i) Add jobs, iff missing
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    nbrOfJobs <- getJobNr(reg);
    if (nbrOfJobs == 0L) {
      verbose && enter(verbose, "Adding jobs to registry");
      more.args <- list(paired=isPaired, indexPrefix=indexPrefix, rgSet=rgSet, params=params, path=outPath, skip=skip, verbose=verbose);
      ids <- batchMap(reg, fun=processOne, getFiles(ds), more.args=more.args);
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
    dW <- 1.00; # Pool every dW seconds
    t0 <- Sys.time();
    tCount <- 0L;
    status <- NULL;
    while (length(findNotTerminated(reg)) > 0L) {
      lastStatus <- status;
      out <- capture.output(status <- showStatus(reg));
      if (identical(status, lastStatus)) {
        verbose && printf(verbose, ".");
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
          verbose && printf(verbose, "[%s]\n", format(dt));
        }
      } else {
        verbose && printf(verbose, "\n");
        verbose && print(verbose, status);
      }
      Sys.sleep(dW);
    } # while(...)
    verbose && exit(verbose);

    verbose && print(verbose, showStatus(reg));

    verbose && exit(verbose);
  } else {
    for (kk in seq_along(ds)) {
      df <- getFile(ds, kk);
      verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", kk, getName(df), length(ds)));

      res <- processOne(df, paired=isPaired, indexPrefix=indexPrefix, rgSet=rgSet, params=params, path=outPath, skip=skip, verbose=verbose);
      verbose && str(verbose, res);

      # Not needed anymore
      df <- res <- NULL;

      verbose && exit(verbose);
    } # for (kk ...)
  }

  res <- getOutputDataSet(this, verbose=less(verbose, 1));

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2013-08-28
# o Now process() outputs distributed status reports on when the status
#   changes.  Inbetween, there is a progress bar.
# 2013-08-26
# o BUG FIX: The internal workaround of process() for Bowtie2Alignment
#   that handled commas in FASTQ filenames deleted the temporary files
#   too soon resulting in garbage alignments.
# o Now process() for Bowtie2Alignment() can utilize BatchJobs.
# 2013-08-24
# o Now Bowtie2Alignment() will do paired-end alignment if the
#   FastqDataSet specifies paired-end reads.
# 2013-08-23
# o BUG FIX: Read Group options ('--rg' and '--rg-id') passed to 'bowtie2'
#   by the Bowtie2Aligment class missed the preceeding '--'.  Also, if
#   the Read Group ID was missing NULL was used - now it is set to 1.
# 2013-07-18
# o Now Bowtie2Alignment handles if there are commas in the pathname of
#   the FASTQ file by using a tempory file link without commas.  This
#   is needed because the bowtie2 executable does not support commas.
# 2013-06-27
# o Now process() for Bowtie2Aligment temporarily decompresses gzipped
#   FASTQ files in case the installed bowtie2 does not support gzip files.
# 2012-11-26
# o BUG FIX: process() for BwaAlignment and Bowtie2Alignment was
#   only align the first sample.
# 2012-10-01
# o Now process() Bowtie2Alignment write SAM read groups, iff given.
# o Now Bowtie2Alignment inherits from AbstractAlignment.
# 2012-09-28
# o Added support for argument 'readGroup' to Bowtie2Alignment().
# 2012-09-27
# o Created from BwaAlignment.R.
############################################################################
