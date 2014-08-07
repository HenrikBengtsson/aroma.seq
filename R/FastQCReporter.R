###########################################################################/**
# @RdocClass FastQCReporter
#
# @title "The FastQCReporter class"
#
# \description{
#  @classhierarchy
#
#  A FastQCReporter takes a @see "FastqDataSet" as input, possibly groups
#  the samples, and generates FastQC [1] reports.
#  How the grouping is done, can be controlled by a parameter.
# }
#
# @synopsis
#
# \arguments{
#  \item{dataSet}{An @see "FastqDataSet".}
#  \item{groupBy}{A @character string or an explicit named @list,
#   specifying which input files should be processed together.}
#  \item{...}{Additional arguments passed to @see "AromaSeqTransform".}
#  \item{.className}{A @character string specifying what class of
#   data sets to accept.}
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
# @author HB
#
# \references{
#  [1] Simon Andrews,
#      FastQC - A quality control tool for high throughput sequence data,
#      March 2014.
#      \url{http://www.bioinformatics.babraham.ac.uk/projects/fastqc/}
# }
#
# \seealso{
#   Internally @see "Rsamtools::mergeBam" is used for merging, sorting,
#   and indexing.
# }
#*/###########################################################################
setConstructorS3("FastQCReporter", function(dataSet=NULL, groupBy=NULL, ..., .className="FastqDataSet") {
  if (!is.null(dataSet)) {
    # Argument 'groupBy':
    if (is.character(groupBy)) {
      groupBy <- match.arg(groupBy, choices=c("name"));
    } else if (is.list(groupBy)) {
      # Validated below
    } else {
      throw("Invalid argument 'groupBy': ", mode(groupBy));
    }
  }

  this <- extend(AromaSeqTransform(dataSet, groupBy=groupBy, ..., .className=.className), c("FastQCReporter", uses("FileGroupsInterface")));

  # Argument 'groupBy':
  if (is.list(groupBy)) {
    validateGroups(this, groups=groupBy);
  }

  this;
})

setMethodS3("getRootPath", "FastQCReporter", function(this, ...) {
  "reports";
}, protected=TRUE)

setMethodS3("getAsteriskTags", "FastQCReporter", function(this, ...) {
  # The 'groupBy' tag
  params <- getParameters(this);
  groupBy <- params$groupBy;
  if (is.character(groupBy)) {
    groupBy <- sprintf("by%s", capitalize(groupBy));
  } else {
    # Generate a short checksum
    groupBy <- getChecksum(groupBy, algo="crc32");
  }

  c("FastQC", groupBy);
}, protected=TRUE)


setMethodS3("getOutputDataSet", "FastQCReporter", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);


  ## Find all existing output data files
  path <- getPath(this);
  fqcs <- FastQCDataFileSet$byPath(path, ...);

  ## Order according to grouped input data set
  names <- getGroupNames(this);
  fqcs <- extract(fqcs, names, onMissing=onMissing);

  fqcs;
}) # getOutputDataSet() for FastQCReporter



setMethodS3("process", "FastQCReporter", function(this, ..., skip=TRUE, force=FALSE, verbose=FALSE) {
  # Requirements
  stopifnot(isCapableOf(aroma.seq, "fastqc"));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Generating FastQC reports");
  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  groups <- getGroups(this);
  verbose && printf(verbose, "Merging into %d groups: %s\n", length(groups), hpaste(names(groups)));
  verbose && str(verbose, head(groups));

  if (force) {
    todo <- seq_along(groups);
  } else {
    todo <- findFilesTodo(this, verbose=less(verbose, 1));
    # Already done?
    if (length(todo) == 0L) {
      verbose && cat(verbose, "Already done. Skipping.");
      res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));
      verbose && exit(verbose);
      return(invisible(res));
    }
  }

  verbose && cat(verbose, "Number of files: ", length(ds));
  verbose && cat(verbose, "Number of groups: ", length(groups));

  path <- getPath(this);
  verbose && cat(verbose, "Output path: ", path);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply FastQC reporting to each group of FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds, IDXS=groups[todo], FUN=function(dfList, ..., force=FALSE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq");

    # Requirements
    stopifnot(isCapableOf(aroma.seq, "fastqc"));

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'dfList':
    stopifnot(is.list(dfList));

    # Argument 'force':
    force <- Arguments$getLogical(force);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    verbose && enter(verbose, "Generating FastQC report");
    verbose && cat(verbose, "Number of FASTQ files to be analyzed: ", length(dfList));

    pathnames <- unlist(sapply(dfList, FUN=getPathname), use.names=FALSE);

    name <- names(dfList)[1L];
    pathD <- file.path(path, name);
    filenameD <- "fastqc_data.txt";
    filenameD <- Arguments$getFilename(filenameD);
    pathnameD <- file.path(pathD, filenameD);
    verbose && cat(verbose, "Output file: ", pathnameD);

    # Already done?
    if (!force && isFile(pathnameD)) {
      verbose && cat(verbose, "Already processed. Skipping.");
      verbose && exit(verbose);
      res <- list(pathnameD=pathnameD);
      return(invisible(res));
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # BEGIN: ATOMIC PROCESSING
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    verbose && enter(verbose, "Atomic processing");
    pathDT <- sprintf("%s.tmp", pathD);
    verbose && cat(verbose, "Temporary output path: ", pathDT);

    # FIXME: Make this method smarter about passing the --casava option;
    # From http://seqanswers.com/forums/showpost.php?p=58152&postcount=200:
    # "Those names don't look like the names generated by Casava. According
    #  to the docs I've got the fastq file names should follow the pattern:
    #
    #  <sample name>_<barcode sequence>_L<lane (0-padded to 3 digits)>_R<read number>_<set number (0-padded to 3 digits>.fastq.gz
    #
    # Which is what FastQC looks for. ..."
    opts <- "--casava";
    verbose && cat(verbose, "FastQC options: ", paste(opts, collapse=", "));
    res <- fastQC(pathnames=pathnames, outPath=pathDT, opts, verbose=less(verbose, 10));
    verbose && cat(verbose, "Output log:");
    verbose && print(verbose, res);

    # Sanity check
    stopifnot(isDirectory(pathDT));


    # FIXME: Do we need to pool here?  Does fastQC() return before output
    # directory is created/available on the file system? /HB 2014-07-18.

    # Scan for output files or output directories (depending on version)

    # (a) FastQC (>= 0.11.1; 2014-06-02)
    #     output is now a self-contained HTML file
    filesT <- list.files(path=pathDT, pattern="_fastqc.html$",
                         ignore.case=TRUE, full.names=FALSE);

    # (b) FastQC (<= 0.10.1; 2012-05-03)
    #     output is a multi-file directory
    dirT <- list.files(path=pathDT, pattern="_fastqc$",
                       ignore.case=TRUE, full.names=FALSE);

    if (length(filesT) > 0L) {
      # FastQC (>= 0.11.1)
      # Nothing more to do...
    } else if (length(dirT) > 0L) {
      # FastQC (<= 0.10.1)
      pathT <- file.path(pathDT, dirT);
      stopifnot(isDirectory(pathT));

      # Sanity check
      pathnameDT <- file.path(pathT, filenameD);
      Arguments$getReadablePathname(pathnameDT);

      # CLEANUP: Remove zip file
      filenameDTZ <- sprintf("%s.zip", dirT);
      pathnameDTZ <- file.path(pathDT, filenameDTZ);
      if (isFile(pathnameDTZ)) {
        file.remove(pathnameDTZ);
      }

      # Move output subdirectory
      file.rename(pathT, pathD);

      # Sanity check
      stopifnot(!isDirectory(pathT));
    } else {
      throw(sprintf("Directories %s contains neither *_fastqc/ directories nor *_fastqc.html files: %s", sQuote(pathDT), paste(sQuote(list.files(path=pathDT, full.names=FALSE)), collapse=", ")));
    }

    # Sanity check
    stopifnot(isDirectory(pathD));

    # CLEANUP
    removeDirectory(pathDT);
    stopifnot(!isDirectory(pathDT));

    verbose && exit(verbose);
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # END: ATOMIC PROCESSING
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    verbose && cat(verbose, "Generated FASTQ report: ", pathD);
    # Sanity checks
    pathD <- Arguments$getReadablePath(pathD);
    pathnameD <- Arguments$getReadablePathname(pathnameD);

    verbose && exit(verbose);

    invisible(list(pathnameD=pathnameD));
  }, force=force, verbose=less(verbose, 1));

  res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 50));

  verbose && exit(verbose);

  invisible(res);
}) # process()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TO DROP
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("validateGroups", "FastQCReporter", function(this, groups, ...) {
  # Input data set
  ds <- getInputDataSet(this);
  nbrOfFiles <- length(ds);

  # Sanity checks
  idxs <- unlist(groups, use.names=FALSE);
  idxs <- Arguments$getIndices(idxs, max=nbrOfFiles);
  if (length(idxs) < nbrOfFiles) {
    throw("One or more input FASTQ files is not part of any group.");
  } else if (length(idxs) > nbrOfFiles) {
    throw("One or more input FASTQ files is part of more than one group.");
  }

  if (is.null(names(groups))) {
    throw("The list of groups does not have names.");
  }

  invisible(groups);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2014-08-07
# o BUG FIX: Now FastQCReporter also works with FastQC (>= 0.11.1).
# 2014-07-18
# o The error message when FastQC fail to generate a '*_fastqc'
#   subdirectory is now more informative.
# 2014-04-13 [HB]
# o CLEANUP: Now FastQCReporter outputs sample directories without the
#   default and auxillary '_fastqc' suffix.
# 2014-03-02 [HB]
# o Created.
############################################################################
