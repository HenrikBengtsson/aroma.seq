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
# @author
#
# \seealso{
#   An object of this class is typically part of an 
#   @see "FastqDataSet".
# }
#*/###########################################################################
setConstructorS3("IlluminaFastqDataFile", function(...) {
  extend(FastqDataFile(...), "IlluminaFastqDataFile");
})



setMethodS3("getSampleName", "IlluminaFastqDataFile", function(this, ...) {
  name <- getFullName(this, ...);
  barcode <- getBarcodeSequence(this);
  pattern <- sprintf("_%s_L.*$", barcode);
  stopifnot(regexpr(pattern, name) != -1);
  name <- gsub(pattern, "", name);
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
    rfq <- readFastq(ff);

    id <- id(rfq)[1L];
    info <- as.character(id);
 
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
      read=gsub(patternB, "\\1", infoB),
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
# 2012-06-29
# o Created.
############################################################################
