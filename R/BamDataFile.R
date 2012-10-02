##########################################################################/**
# @RdocClass BamDataFile
#
# @title "The abstract BamDataFile class"
#
# \description{
#  @classhierarchy
#
#  A BamDataFile object represents a BAM file.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#
# \references{
#  [1] The SAM Format Specification Working Group,
#      \emph{The SAM Format Specification}, Sept 7, 2011.\cr
# }
#
# \seealso{
#   An object of this class is typically part of an 
#   @see "BamDataSet".
# }
#*/###########################################################################
setConstructorS3("BamDataFile", function(...) {
  extend(GenericDataFile(...), "BamDataFile");
})


setMethodS3("as.character", "BamDataFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  s <- c(s, sprintf("Has index file (*.bai): %s", hasIndex(this)));
  s <- c(s, sprintf("Is sorted: %s", isSorted(this)));
  if (hasIndex(this)) {
    n <- nbrOfTargets(this);
    s <- c(s, sprintf("Number of targets: %s", n));
    len <- getTotalTargetLength(this);
    s <- c(s, sprintf("Total target length: %.3gMb (%.0f bases)", len/1e9, len));
    names <- getTargetNames(this);
    s <- c(s, sprintf("Targets: [%d] %s", n, hpaste(names)));

    counts <- getReadCounts(this);
    total <- sum(counts, na.rm=TRUE);
    s <- c(s, sprintf("Number of mapped reads: %d (%.1f%%) out of %d", counts[["mapped"]], 100*counts[["mapped"]]/total, total));
    s <- c(s, sprintf("Number of unmapped reads: %d (%.1f%%) out of %d", counts[["unmapped"]], 100*counts[["unmapped"]]/total, total));
  }

  class(s) <- class;
  s;
})

setMethodS3("buildIndex", "BamDataFile", function(this, ..., skip=!overwrite, overwrite=FALSE) {
  pathname <- getPathname(this);
  pathnameBAI <- sprintf("%s.bai", pathname);

  if (hasIndex(this)) {
    if (skip) {
      return(invisible(getIndexFile(this)));
    }
    if (!overwrite) {
      throw("Cannot build index file (*.bai). File already exists: ", pathnameBAI);
    }
  }

  pathnameT <- Rsamtools::indexBam(pathname);

  getIndexFile(this);
})


setMethodS3("getIndexFile", "BamDataFile", function(this, ...) {
  pathname <- getPathname(this);
  pathnameBAI <- sprintf("%s.bai", pathname);
  if (!isFile(pathnameBAI)) {
    return(NULL);
  }
  BamIndexDataFile(pathnameBAI);
})


setMethodS3("hasIndex", "BamDataFile", function(this, ...) {
  !is.null(getIndexFile(this));
})



setMethodS3("getIndexStats", "BamDataFile", function(this, ..., force=FALSE) {
  stats <- this$.idxStats;

  if (force || is.null(stats)) {
    pathname <- getPathname(this);
    if (!hasIndex(this)) {
      throw("Cannot get index statistics. Index does not exists: ", pathname);
    }
    bfr <- systemSamtools("idxstats", pathname, stdout=TRUE, ...);
    bfr <- strsplit(bfr, split="\t", fixed=TRUE);
    seqName <- sapply(bfr, FUN=.subset, 1L);
    seqLength <- as.integer(sapply(bfr, FUN=.subset, 2L));
    countMapped <- as.integer(sapply(bfr, FUN=.subset, 3L));
    countUnmapped <- as.integer(sapply(bfr, FUN=.subset, 4L));

    # NOTE: samtools idxstats can return ridicolously(!) large read
    # counts.  Those will become NAs in as.integer() above.  Let's
    # assume they are errors. /HB 2010-10-02

    stats <- data.frame(length=seqLength, mapped=countMapped, unmapped=countUnmapped);
    rownames(stats) <- seqName;
    this$.idxStats <- stats;
  }

  stats;
})


setMethodS3("getReadCounts", "BamDataFile", function(this, ...) {
  stats <- getIndexStats(this, ...);
  counts <- colSums(stats[,c("mapped", "unmapped")], na.rm=TRUE);
  counts;
})


setMethodS3("nbrOfReads", "BamDataFile", function(this, ...) {
  counts <- getReadCounts(this, ...);
  sum(counts, na.rm=TRUE);
})

setMethodS3("nbrOfMappedReads", "BamDataFile", function(this, ...) {
  counts <- getReadCounts(this, ...);
  counts["mapped"];
})

setMethodS3("nbrOfUnmappedReads", "BamDataFile", function(this, ...) {
  counts <- getReadCounts(this, ...);
  counts["unmapped"];
})


# \details{
#   BAM headers typically contain an \code{"@HD VN:1.0 SO:<value>"} entry,
#   where \code{<value>} indicates whether the aligned reads are sorted 
#   or not.  Unfortunately, this entry is neither enforced nor has it to
#   be correct [1,2].
#
#   Instead, we consider a BAM file to be sorted if and only if it has
#   an index file.  The rationale is that it is not possible to index
#   a BAM file unless it is sorted first.
# }
#
# \references{
#   [1] Question: is my BAM file sorted?, Biostar, 2011,
#   \url{http://www.biostars.org/post/show/5256/is-my-bam-file-sorted/}\cr
#   [2] Asking for suggestiona on samtools bug fixing, SEQanswers, 2010,
#   \url{http://seqanswers.com/forums/showthread.php?t=3739}\cr
# }
setMethodS3("isSorted", "BamDataFile", function(this, ...) {
  isTRUE(hasIndex(this));
})


# Argument '...' must be 2nd to match the generic base::sort() function.
setMethodS3("sort", "BamDataFile", function(x, ..., force=FALSE) {
  # To please R CMD check
  this <- x;

  # Nothing todo?
  if (!force && isSorted(this)) {
    return(this);
  }

  throw("Not yet implemented!");
})


setMethodS3("nbrOfSeqs", "BamDataFile", function(this, ...) {
  nbrOfTargets(this);
})

setMethodS3("getTargets", "BamDataFile", function(this, ...) {
  hdr <- getHeader(this);
  targets <- hdr$targets;
  targets;  
})

setMethodS3("nbrOfTargets", "BamDataFile", function(this, ...) {
  length(getTargets(this));
})

setMethodS3("getTargetNames", "BamDataFile", function(this, ...) {
  names(getTargets(this));
})

setMethodS3("getTargetLengths", "BamDataFile", function(this, ...) {
  getTargets(this);
})

setMethodS3("getTotalTargetLength", "BamDataFile", function(this, ...) {
  sum(as.numeric(getTargets(this)));
})

setMethodS3("getHeader", "BamDataFile", function(this, force=FALSE, ...) {
  header <- this$.header;
  if (force || is.null(header)) {
    header <- readHeader(this, ...);
    this$.header <- header;
  }
  header;
})

setMethodS3("readHeader", "BamDataFile", function(this, ...) {
  pathname <- getPathname(this);
  if (!hasIndex(this)) {
    throw("Cannot read header. Index file (*.bai) is missing: ", pathname);
  }
  bf <- Rsamtools::BamFile(pathname);
  hdr <- Rsamtools::scanBamHeader(bf);
  hdr;
}, private=TRUE)


# \section{Read Groups (RG) specification}{
#   The RG fields populated are ID, SM, PL and LB, where
#   ID=identifier, SM=sample, PL=platform/technology, and
#   LB=library (DNA library prep identifier).
#   ID: Each RG must have a unique ID.
#   SM: Sample. Use pool name where a pool is being sequenced.
#   PM: Platform/technology used to produce the reads. Valid values: 
#       CAPILLARY, HELICOS, ILLUMINA, IONTORRENT, LS454, PACBIO, and SOLID.
#   LB: Library.
#
#   Note that tools such as GATK (Broad Institute) requires RG:s.
# }
#
# \references{
#  [1] The SAM Format Specification Working Group,
#      \emph{The SAM Format Specification}, Sept 7, 2011.\cr
#  [2] "What is '@RG ID' versus '@RG SM'"
#      \url{http://seqanswers.com/forums/showthread.php?t=9784}\cr
# }
setMethodS3("getReadGroups", "BamDataFile", function(this, ...) {
  hdr <- getHeader(this, ...);
  SamReadGroup$byScanBamHeader(hdr, ...);
})

setMethodS3("getReadGroup", "BamDataFile", function(this, ...) {
  rgList <- getReadGroups(this, ...);
  if (length(rgList) >= 1L) {
    rgList <- rgList[[1L]];
  } else {
    rgList <- SamReadGroup();
  }
  rgList;
})


# @RdocMethod replaceAllReadGroups
# @title "Writes a new BAM file with all existing read groups replaced by one new read group"
#
# \description{
#   @get "title".
# }
#
# \arguments{
#  \item{sample}{Specifies the \code{SM} read group.}
#  \item{library}{Specifies the \code{LB} read group.}
#  \item{platform}{Specifies the \code{PL} read group.}
#  \item{platformUnit}{Specifies the \code{PU} read group.}
# }
setMethodS3("replaceAllReadGroups", "BamDataFile", function(this, rg="*", ..., validate=TRUE, skip=!overwrite, overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'rg':
  if (is.character(rg) && rg == "*") {
    rg <- SamReadGroup();
    keys <- names(asSamList(rg, drop=FALSE));
    for (key in keys) rg[[key]] <- "*";
  } else {
    rg <- Arguments$getInstanceOf(rg, "SamReadGroup");
  }


  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating Read Groups");
 

  pathname <- getPathname(this);

  # Output pathname
  pathnameD <- gsub("[.]bam$", ",RG.bam", pathname);
  pathnameD <- Arguments$getWritablePathname(pathnameD, mustNotExist=(!skip && !overwrite));

  if (isFile(pathnameD)) {
    if (skip) {
      verbose && cat(verbose, "Already exists. Skipping.");
      bf <- newInstance(this, pathnameD);
      buildIndex(bf, skip=TRUE, verbose=less(verbose, 10));
      verbose && exit(verbose);
      return(bf);
    }
    if (!overwrite) {
      throw("Cannot replace SAM read groups. File already exists: ", pathnameD);
    }
  }


  verbose && cat(verbose, "Arguments:");
  verbose && print(verbose, rg);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Use default values from existing read groups of the input file?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rgList <- asSamList(rg);
  keep <- sapply(rgList, FUN=function(x) identical(x, "*"));
  keys <- names(rgList)[keep];
  if (length(keys) > 0L) {
    rgDefault <- getReadGroups(this);
    rgDefault <- Reduce(merge, rgDefault);
    # Sanity check
    stopifnot(inherits(rgDefault, "SamReadGroup"));

    for (key in keys) {
      rg[[key]] <- rgDefault[[key]];
    }
  } # for (key ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  validate(rg);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert mandatory fields 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rgList <- asSamList(rg);
  mandatory <- c("SM", "LB", "PL", "PU");
  missing <- setdiff(mandatory, names(rgList));
  if (length(missing)) {
    throw("Cannot write read groups. Mandatory fields are missing: ", hpaste(missing));
  }

  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, rgList);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write new BAM file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list("AddOrReplaceReadGroups", I=pathname, O=pathnameD);
  args <- c(args, asString(rg, fmtstr="RG%s=%s"));
  verbose && str(verbose, args);
  args$verbose <- less(verbose, 10);

  res <- do.call(systemPicard, args);

  bf <- newInstance(this, pathnameD);
  buildIndex(bf, overwrite=TRUE, verbose=less(verbose, 10));

  verbose && exit(verbose);

  bf;
}) # replaceAllReadGroups()


############################################################################
# HISTORY:
# 2012-10-02
# o Added getIndexStats() for BamDataFile.
# o Now buildIndex() for BamDataFile returns a BamIndexDataFile.
# o Added getIndexFile() for BamDataFile.
# 2012-06-28
# o Added getReadGroups() and replaceAllReadGroups() and buildIndex().
# o Created.
############################################################################
