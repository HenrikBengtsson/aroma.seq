###########################################################################/**
# @RdocClass IlluminaFastqDataFile
#
# @title "The abstract IlluminaFastqDataFile class"
#
# \description{
#  @classhierarchy
#
#  A IlluminaFastqDataFile object represents a FASTQ data file.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "FastqDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# \seealso{
#   An object of this class is typically part of an
#   @see "FastqDataSet".
# }
#*/###########################################################################
setConstructorS3("IlluminaFastqDataFile", function(...) {
  extend(FastqDataFile(...), "IlluminaFastqDataFile");
})



setMethodS3("getFileVersion", "IlluminaFastqDataFile", function(this, ...) {
  name <- getFullName(this, ...);
  patterns <- c("Casava_v1.4"="^[^_]+_[ACGTN]+_L[0-9]+_R[0-9]");
  for (key in names(patterns)) {
    pattern <- patterns[key];
    if (regexpr(pattern, name) != -1) {
      return(key);
    }
  }
  as.character(NA);
})


setMethodS3("getSampleName", "IlluminaFastqDataFile", function(this, ...) {
  name <- getFullName(this, ...);

  ver <- getFileVersion(this);
  if (is.na(ver)) ver <- "<gzipped; unknown>";
  if (ver == "Casava_v1.4") {
    barcode <- getBarcodeSequence(this);
    # AD HOC patch for observing ATGNCA when expected ATGTCA. /HB 2012-10-01
    barcode <- gsub("N", ".", barcode, fixed=TRUE);
    pattern <- sprintf("_%s_L[0-9]+_R[0-9](_[0-9]+)$", barcode);
    if (regexpr(pattern, name) == -1L) {
      throw(sprintf("The fullname (%s) of the %s with version %s does not match the expected pattern (%s): %s", sQuote(name), class(this)[1L], sQuote(ver), sQuote(pattern), getPathname(this)));
    }
    name <- gsub(pattern, "", name);
  } else {
    warning("Unknown Illumina FASTQ file version. Using fullname as sample name: ", name);
  }

  name;
})

setMethodS3("getPlatform", "IlluminaFastqDataFile", function(this, ...) {
  "Illumina";
})

setMethodS3("getLane", "IlluminaFastqDataFile", function(this, ...) {
  info <- getFirstSequenceInfo(this);
  info$laneIdx;
})

setMethodS3("getInstrumentId", "IlluminaFastqDataFile", function(this, ...) {
  info <- getFirstSequenceInfo(this);
  info$instrumentId;
})

setMethodS3("getFlowcellId", "IlluminaFastqDataFile", function(this, ...) {
  info <- getFirstSequenceInfo(this);
  info$flowcellId;
})

setMethodS3("getBarcodeSequence", "IlluminaFastqDataFile", function(this, ...) {
  info <- getFirstSequenceInfo(this);
  info$indexSequence;
})

setMethodS3("getReadDirection", "IlluminaFastqDataFile", function(this, ...) {
  info <- getFirstSequenceInfo(this);
  info$read;
})

setMethodS3("getPlatformUnit", "IlluminaFastqDataFile", function(this, ...) {
#    PU: the "platform unit" - a unique identifier which tells you what
#        run/experiment created the data.  For Illumina, please follow this
#        convention: Illumina flowcell barcode suffixed with a period and
#        the lane number (and further suffixed with period followed by
#        sample member name for pooled runs). If referencing an existing
#        already archived run, then please use the run alias in the SRA.
  parts <- c(getFlowcellId(this), getLane(this), getSampleName(this));
  paste(parts, collapse=".");
})

setMethodS3("getFirstSequenceInfo", "IlluminaFastqDataFile", function(this, force=FALSE, ...) {
  require("ShortRead") || throw("Package not loaded: ShortRead");

  info <- this$.info;

  if (force || is.null(info)) {
    pathnameFQ <- getPathname(this);
    ff <- FastqFile(pathnameFQ);
    on.exit(close(ff));

    # Really inefficient way to find the first sequence information.
    # /HB 2013-11-19
    rfq <- readFastq(ff);

    id <- id(rfq)[1L];
    info <- as.character(id);
    rfq <- NULL; # Not needed anymore

    patternA <- "^([^:]+):([0-9]+):([^:]+):([0-9]+):([0-9]+):([0-9]+):([0-9]+)";
    patternB <- " ([^:]+):([^:]+):([0-9]+):([^:]+)$";
    pattern <- sprintf("%s%s", patternA, patternB);
    stopifnot(regexpr(pattern, info) != -1);

    infoA <- gsub(patternB, "", info);
    infoB <- gsub(patternA, "", info);

    info <- list(
      instrumentId=gsub(patternA, "\\1", infoA),
      runIdx=as.integer(gsub(patternA, "\\2", infoA)),
      flowcellId=gsub(patternA, "\\3", infoA),
      laneIdx=as.integer(gsub(patternA, "\\4", infoA)),
      tileIdx=as.integer(gsub(patternA, "\\5", infoA)),
      x=as.integer(gsub(patternA, "\\6", infoA)),
      y=as.integer(gsub(patternA, "\\7", infoA)),
      read=as.integer(gsub(patternB, "\\1", infoB)),
      isFiltered=gsub(patternB, "\\2", infoB),
      controlNumber=as.integer(gsub(patternB, "\\3", infoB)),
      indexSequence=gsub(patternB, "\\4", infoB)
    );

    this$.info <- info;
  }

  info;
}, protected=TRUE)



setMethodS3("getDefaultSamReadGroup", "IlluminaFastqDataFile", function(this, ...) {
  # SM: Sample
  # PL: Platform unit
  # PU: Platform
  SM <- getSampleName(this);
  PL <- getPlatform(this);
  PU <- getPlatformUnit(this);
  SamReadGroup(SM=SM, PL=PL, PU=PU);
})



############################################################################
# HISTORY:
# 2013-11-19
# o Now getSampleName() for IlluminaFastqDataFile gives a more informative
#   error message if the fullname does not match the expected pattern.
# 2012-10-16
# o BUG FIX: getSampleName() did not handle trailing indices, e.g. _002.
# 2012-10-02
# o Added getFileVersion().
# 2012-06-29
# o Created.
############################################################################
