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


  bowtie2 <- function(pathnameFQ, indexPrefix, pathnameSAM, ..., gzAllowed=NA, verbose=FALSE) {
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
  params <- getParameters(this);
  verbose && cat(verbose, "Additional bowtie2 arguments:");
  verbose && str(verbose, params);

  verbose && cat(verbose, "Number of files: ", length(ds));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds, FUN=function(df, paired=FALSE, indexPrefix, rgSet, params, path, ...., skip=TRUE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'df':
    df <- Arguments$getInstanceOf(df, "FastqDataFile");

    # Argument 'paired':
    paired <- Arguments$getLogical(paired);

    # Argument 'indexPrefix':
    indexPrefix <- Arguments$getCharacter(indexPrefix);

    # Argument 'skip':
    skip <- Arguments$getLogical(skip);

    # Argument 'path':
    path <- Arguments$getWritablePath(path);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    verbose && enter(verbose, "Bowtie2 alignment of one sample");

    # Paired-end?
    if (paired) {
      pathnameFQ <- sapply(list(R1=df, R2=getMateFile(df)), FUN=getPathname);
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
    } # if (done)

    verbose && exit(verbose);

    invisible(list(pathnameFQ=pathnameFQ, pathnameSAM=pathnameSAM, pathnameBAM=pathnameBAM));
  }, paired=isPaired(this), indexPrefix=indexPrefix, rgSet=rgSet, params=params, path=getPath(this), skip=skip, verbose=verbose) # dsApply()

  res <- getOutputDataSet(this, verbose=less(verbose, 1));

  verbose && exit(verbose);

  invisible(res);
})


############################################################################
# HISTORY:
# 2013-08-31
# o Now process() for Bowtie2Alignment utilizes dsApply().
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
