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

    rgArgs <- asString(rg, fmtstr="%s:%s");
    rgArgs <- rgArgs[regexpr("^ID:", rgArgs) == -1];

    # Don't forget to put within quotation marks
    rgArgs <- sprintf("\"%s\"", rgArgs);

    rgArgs <- as.list(rgArgs);
    names(rgArgs) <- rep("rg", times=length(rgArgs));

    rgArgs <- c(list("rg-id"=asSamList(rg)$ID), rgArgs);

    rgArgs;
  } # asBowtie2Parameters()


  bowtie2 <- function(pathnameFQ, indexPrefix, pathnameSAM, ..., verbose=FALSE) {
    pathnameFQ <- Arguments$getReadablePathname(pathnameFQ);
    indexPrefix <- Arguments$getCharacter(indexPrefix);
    indexPath <- Arguments$getReadablePath(getParent(indexPrefix));
    pathnameSAM <- Arguments$getWritablePathname(pathnameSAM);
    systemBowtie2(args=list("-x"=indexPrefix, "-U"=pathnameFQ, "-S"=pathnameSAM, ...), verbose=verbose);
  } # bowtie2()



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

  outPath <- getPath(this);
  for (kk in seq_along(ds)) {
    df <- getFile(ds, kk);
    name <- getName(df);
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d",
                                                    kk, name, length(ds)));

    pathnameFQ <- getPathname(df);
    verbose && cat(verbose, "FASTQ pathname: ", pathnameFQ);

    # The SAM and BAM files to be generated
    fullname <- getFullName(df);
    filename <- sprintf("%s.sam", fullname);
    pathnameSAM <- Arguments$getWritablePathname(filename, path=outPath);
    verbose && cat(verbose, "SAM pathname: ", pathnameSAM);
    filename <- sprintf("%s.bam", fullname);
    pathnameBAM <- Arguments$getWritablePathname(filename, path=outPath);
    verbose && cat(verbose, "BAM pathname: ", pathnameBAM);

    # Nothing to do?
    if (skip && isFile(pathnameBAM)) {
      verbose && cat(verbose, "Already aligned. Skipping");
      verbose && exit(verbose);
      next;
    }

    verbose && print(verbose, df);

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # (a) Generate SAM file
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!isFile(pathnameSAM)) {
      if (isGzipped(df)) {
        if (is.na(gzAllowed)) {
          gzAllowed <- queryBowtie2("support:fastq.gz");
        }
        if (!gzAllowed) {
          why <- attr(gzAllowed, "why");
          throw(sprintf("Cannot align reads in '%s': %s", getPathname(df), why));
        }
      }

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

      # In case bowtie2 generates empty SAM files
      # /HB 2012-10-01 (still to be observed)
      if (isFile(pathnameSAM)) {
        if (file.info(pathnameSAM)$size == 0L) {
          verbose && cat(verbose, "Removing empty SAM file falsely created by Bowtie2: ", pathnameSAM);
          file.remove(pathnameSAM);
        }
      }

    }

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

    verbose && exit(verbose);
  } # for (kk ...)

  res <- getOutputDataSet(this, verbose=less(verbose, 1));

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
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
