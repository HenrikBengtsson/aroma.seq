###########################################################################/**
# @RdocClass FastqDataFile
#
# @title "The abstract FastqDataFile class"
#
# \description{
#  @classhierarchy
#
#  A FastqDataFile object represents a FASTQ data file.
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
# @author "HB"
#
# \references{
#   Wikipedia, FASTQ format,
#   \url{http://en.wikipedia.org/wiki/FASTQ_format}.\cr
# }
#
# \seealso{
#   An object of this class is typically part of an
#   @see "FastqDataSet".
# }
#*/###########################################################################
setConstructorS3("FastqDataFile", function(...) {
  extend(GenericDataFile(...), "FastqDataFile");
})


setMethodS3("as.character", "FastqDataFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  class <- class(s);

  n <- nbrOfSeqs(this);
  s <- c(s, sprintf("Number of sequences: %s", n));
  s <- c(s, sprintf("Common width of sequences: %d", getCommonSeqWidth(this)));

  class(s) <- class;
  s;
}, protected=TRUE)


setMethodS3("getDefaultFullName", "FastqDataFile", function(this, ...) {
  name <- NextMethod("getDefaultFullName");
  name <- gsub("[.](fastq|fq)$", "", name, ignore.case=TRUE);
  name;
}, protected=TRUE)


setMethodS3("nbrOfSeqs", "FastqDataFile", function(this, ...) {
  geo <- getGeometry(this, ...);
  geo[1L];
})


setMethodS3("getCommonSeqWidth", "FastqDataFile", function(this, ...) {
  geo <- getGeometry(this, ...);
  geo[2L];
})


setMethodS3("getGeometry", "FastqDataFile", function(this, force=FALSE, ...) {
  geometry <- this$.geometry;
  if (force || is.null(geometry)) {
    geometry <- readGeometry(this, ...);
    if (!anyMissing(geometry)) {
      this$.geometry <- geometry;
    }
  }
  geometry;
})

setMethodS3("readGeometry", "FastqDataFile", function(this, ...) {
  # TO DO: Support gzipped files. /HB 2013-06-20
  if (isGzipped(this)) {
    return(c(NA_integer_, NA_integer_));
  }
  pathname <- getPathname(this);
  geo <- Biostrings::fastq.geometry(pathname);
  geo;
}, private=TRUE)


setMethodS3("getDefaultSamReadGroup", "FastqDataFile", function(this, ...) {
  SamReadGroup();
})

setMethodS3("setSamReadGroup", "FastqDataFile", function(this, rg, ...) {
  # Argument 'rg':
  if (!is.null(rg)) {
    rg <- Arguments$getInstanceOf(rg, "SamReadGroup");
  }
  this$.rg <- rg;
  invisible(this);
})

setMethodS3("getSamReadGroup", "FastqDataFile", function(this, ...) {
  rg <- this$.rg;
  if (is.null(rg)) {
    rg <- getDefaultSamReadGroup(this, ...);
  }
  rg;
})


############################################################################
# HISTORY:
# 2013-06-25
# o Added getDefaultFullName() for FastqDataFile so <fullname>.fastq.gz
#   is properly handled.  Should ideally handled by R.filesets.
# 2013-06-20
# o Now readGeometry() and getGeometry() returns c(NA,NA) for
#   gzipped files. This is better than an error.
# 2012-06-28
# o Created.
############################################################################
