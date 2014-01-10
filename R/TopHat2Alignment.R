###########################################################################/**
# @RdocClass TopHat2Alignment
#
# @title "The TopHat2Alignment class"
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
#  \item{groupBy}{A @character string or an explicit named @list,
#   specifying which input files should be processed together.}
#  \item{indexSet}{An @see "Bowtie2IndexSet".}
#  \item{transcripts}{Gene model (transcriptome) GTF/GFF3 file.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Supported operating systems}{
#   This method is available on Linux and OSX [1].
# }
#
# @author "TT"
#
# \references{
#  [1] TopHat, University of Maryland, 2013.
#      \url{http://http://tophat.cbcb.umd.edu/} \cr
#  [2] Trapnell et al. \emph{Differential gene and transcript expression
#      analysis of RNA-seq experiments with TopHat and Cufflinks}.
#      Nat Protoc, 2012.\cr
# }
#*/###########################################################################
setConstructorS3("TopHat2Alignment", function(..., groupBy=NULL, indexSet=NULL, transcripts=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'groupBy':
  if (is.null(groupBy)) {
  } else if (is.character(groupBy)) {
    groupBy <- match.arg(groupBy, choices=c("name"));
  } else if (is.list(groupBy)) {
    # Validated below
  } else {
    throw("Invalid argument 'groupBy': ", mode(groupBy));
  }

  # Argument 'indexSet':
  if (!is.null(indexSet)) {
    indexSet <- Arguments$getInstanceOf(indexSet, "Bowtie2IndexSet");
  }

  # Argument 'transcripts':
  if (!is.null(transcripts)) {
    transcripts <- Arguments$getInstanceOf(transcripts, "GenericDataFile");
  }

  # Arguments '...':
  args <- list(...);

  this <- extend(AbstractAlignment(..., indexSet=indexSet), "TopHat2Alignment",
    transcripts = transcripts,
    groupBy = groupBy
  );

  # Argument 'groupBy':
  if (is.list(groupBy)) {
    validateGroups(this, groups=groupBy);
  }

  this;
})


setMethodS3("getRootPath", "TopHat2Alignment", function(this, ...) {
  "tophat2Data";
}, protected=TRUE)



setMethodS3("validateGroups", "TopHat2Alignment", function(this, groups, ...) {
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


setMethodS3("getGroups", "TopHat2Alignment", function(this, ...) {
  ds <- getInputDataSet(this);

  groups <- this$groupBy;
  if (is.null(groups)) {
    groups <- as.list(seq_along(ds));
    names(groups) <- getFullNames(ds);
  } else if (is.character(groups)) {
    if (groups == "name") {
      names <- getNames(ds);
      unames <- unique(names);
      idxs <- match(names, unames);
      names(idxs) <- getFullNames(ds);
      groups <- tapply(idxs, INDEX=idxs, FUN=list);
      names(groups) <- unames;
    }
  }
  # Sanity check
  stopifnot(is.list(groups));

  # Range check
  max <- length(ds);
  groups <- lapply(groups, FUN=Arguments$getIndices, max=max);

  # Names
  if (is.null(names(groups))) {
    names(groups) <- sprintf("Group_%d", seq_along(groups));
  }

  groups;
}, protected=TRUE) # getGroups()


setMethodS3("nbrOfGroups", "TopHat2Alignment", function(this, ...) {
  length(getGroups(this));
})


setMethodS3("getOutputDataSet", "TopHat2Alignment", function(this, onMissing=c("drop", "NA", "error"), ...) {
  onMissing <- match.arg(onMissing);
  path <- getPath(this);
  bams <- BamDataSet$byPath(path, ...);
  groups  <- getGroups(this);
  fullnames <- names(groups);
  bams <- extract(bams, fullnames, onMissing=onMissing);
  bams
})


setMethodS3("getSampleNames", "TopHat2Alignment", function(this, ...) {
  ds <- getInputDataSet(this);
  sampleNames <- sub("_(1|R1)$", "", getFullNames(ds));
  sampleNames;
})

setMethodS3("getExpectedOutputPaths", "TopHat2Alignment", function(this, ...) {
  # Find all available output directories
  path <- getPath(this);
  sampleNames <- getSampleNames(this);
  paths <- file.path(path, sampleNames);
  paths;
}, protected=TRUE)


setMethodS3("process", "TopHat2Alignment", function(this, ..., skip=TRUE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Test for non-compatible bowtie2 and tophat2 versions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  stopifnot(isCapableOf(aroma.seq, "tophat2"));
  stopifnot(isCapableOf(aroma.seq, "bowtie2"));
  verT <- attr(findTopHat2(), "version");
  verB <- attr(findBowtie2(), "version");
  if (!is.null(verT) && !is.null(verB)) {
    bad <- (verT == "2.0.3" && verB == "2.1.0");
    if (bad) {
      throw(sprintf("Detected incompatible software installations. TopHat2 v%s is known to not work with Bowtie2 v%s.", verT, verB))
    }
  }
  verT <- verB <- NULL;


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


  verbose && enter(verbose, "TopHat2 alignment");
  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get groups of items to be processed at the same time
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Grouping input data set");
  groups <- getGroups(this);
  verbose && printf(verbose, "Merging into %d groups: %s\n", length(groups), hpaste(names(groups)));
  verbose && str(verbose, head(groups));
  verbose && cat(verbose, "Number of items per groups:");
  ns <- sapply(groups, FUN=length);
  t <- table(ns);
  names(t) <- sprintf("n=%s", names(t));
  verbose && print(verbose, t);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify groups to be processed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (force) {
    todo <- seq_along(groups);
  } else {
    bams <- getOutputDataSet(this, onMissing="NA", verbose=less(verbose, 1));
    todo <- which(!sapply(bams, FUN=isFile));
  }
  verbose && cat(verbose, "Number of groups to process: ", length(todo));

  # Already done?
  if (!force && length(todo) == 0L) {
    verbose && cat(verbose, "Already processed.");
    verbose && print(verbose, bams);
    verbose && exit(verbose);
    return(bams);
  }

  isPaired <- isPaired(ds);
  verbose && cat(verbose, "Paired-end analysis: ", isPaired);

  is <- getIndexSet(this);
  verbose && cat(verbose, "Aligning using index set:");
  verbose && print(verbose, is);

  outPath <- getPath(this);
  verbose && cat(verbose, "Output directory: ", outPath);

  transcripts <- this$transcripts;
  verbose && cat(verbose, "Using transcripts:");
  verbose && print(verbose, transcripts);

  verbose && cat(verbose, "Number of files: ", length(ds));
  verbose && cat(verbose, "Number of groups: ", length(groups));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds, IDXS=groups[todo], DROP=FALSE, FUN=function(dfListR1, isPaired=FALSE, indexSet, transcripts=NULL, outPath, ...., skip=TRUE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq, Rsamtools");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'dfListR1':
    dfListR1 <- Arguments$getInstanceOf(dfListR1, "list");
    dfListR1 <- Arguments$getVector(dfListR1, length=c(1,Inf));
    dfListR1 <- lapply(dfListR1, FUN=Arguments$getInstanceOf, "FastqDataFile");

    # Argument 'isPaired':
    isPaired <- Arguments$getLogical(isPaired);

    # Argument 'indexSet':
    indexSet <- Arguments$getInstanceOf(indexSet, "Bowtie2IndexSet");

    # Argument 'skip':
    skip <- Arguments$getLogical(skip);

    # Argument 'outPath':
    outPath <- Arguments$getWritablePath(outPath);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    sampleName <- sub("_(1|R1)$", "", getFullName(dfListR1[[1L]]));
    verbose && enter(verbose, "Sample name ", sQuote(sampleName));

    gtf <- NULL;
    if (!is.null(transcripts)) gtf <- getPathname(transcripts);

    reads1 <- sapply(dfListR1, FUN=getPathname);
    verbose && printf(verbose, "R1 FASTQ files: [%d] %s\n", length(reads1), hpaste(sQuote(reads1)));

    # Final sample-specific output directory
    outPathS <- file.path(outPath, sampleName);
    args <- list(
      bowtieRefIndexPrefix=getIndexPrefix(indexSet),
      reads1=reads1,
      reads2=NULL,
      outPath=outPathS,
      gtf=gtf
    );

    if (isPaired) {
      dfListR2 <- lapply(dfListR1, FUN=getMateFile);
      reads2 <- sapply(dfListR2, FUN=getPathname);
      verbose && printf(verbose, "R2 FASTQ files: [%d] %s\n", length(reads2), hpaste(sQuote(reads2)));
      args$reads2 <- reads2;
    }


    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # BEGIN: ATOMIC OUTPUT
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Write to temporary output directory
    args$outPath <- sprintf("%s.tmp", args$outPath);
    verbose && cat(verbose, "Temporary output directory: ", args$outPath);

    # (a) Align reads using TopHat2
    verbose && cat(verbose, "Arguments passed to TopHat:");
    verbose && str(verbose, args);
    args$verbose <- less(verbose, 1);
    res <- do.call(tophat2, args=args);

    # (b) Generates BAM index file (assuming the BAM file is sorted)
    pathnameBAM <- file.path(args$outPath, "accepted_hits.bam");
    verbose && cat(verbose, "BAM file: ", pathnameBAM);
    pathnameBAI <- indexBam(pathnameBAM);
    verbose && cat(verbose, "BAM index file: ", pathnameBAI);

    # Rename from temporary to final directory
    file.rename(args$outPath, outPathS);
    verbose && cat(verbose, "Final output directory: ", outPathS);
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # END: ATOMIC OUTPUT
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

    verbose && exit(verbose);

    invisible(list(res=res));
  }, isPaired=isPaired, indexSet=is, transcripts=transcripts, outPath=getPath(this), skip=skip, verbose=verbose) # dsApply()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bams <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));
  verbose && print(verbose, bams);

  # Sanity check
  stopifnot(all(sapply(bams, FUN=isFile)));

  verbose && exit(verbose);

  bams;
})


############################################################################
# HISTORY:
# 2014-01-09 [HB]
# o Added support for argument 'groupBy' to TopHat2Alignment.
# 2013-11-21 [HB]
# o Now process() for TopHat2Alignment generates an index file for
#   accepted_hits.bam.  This will also assert the assumption that
#   TopHat2 outputs sorted BAM files, because if not indexing will
#   fail and an error will be generated.
# 2013-11-16 [HB]
# o CLEANUP: Dropped several methods now taken care of by super class.
# 2013-11-01 [HB]
# o SPEEDUP: Parallized process() for TopHat2Alignment.
# o Now process() for TopHat2Alignment skips already processed samples.
# o Now process() for TopHat2Alignment should also work for single-end
#   reads as well as paired-end reads.
# 2013-10-31 [HB]
# o Now utilizing tophat2().
# o CLEANUP: Dropped non-used argument 'outputDir' from constructor.
# o CLEANUP: Dropping stray cut'n'paste code.
# 2013-10-30 [HB]
# o ROBUSTNESS: Now process() for TopHat2Alignment checks for known
#   incompatible versions of bowtie2 and tophat2.
# 2013-08-12 [TT]
# o Created from Bowtie2Alignment.R.
############################################################################
