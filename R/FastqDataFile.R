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
#   An object of this class is typically part of a @see "FastqDataSet".
# }
#*/###########################################################################
setConstructorS3("FastqDataFile", function(...) {
  extend(GenericDataFile(...), "FastqDataFile");
})


setMethodS3("as.character", "FastqDataFile", function(x, ...) {
  this <- x;
  s <- NextMethod("as.character");
  n <- nbrOfSeqs(this);
  s <- c(s, sprintf("Number of sequences: %s", n));
  s <- c(s, sprintf("Common width of sequences: %d", getCommonSeqWidth(this)));
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
    if (!Biobase::anyMissing(geometry)) {
      this$.geometry <- geometry;
    }
  }
  geometry;
})


setMethodS3("readGeometry", "FastqDataFile", function(this, ...) {
  naValue <- c(NA_integer_, NA_integer_);

  # Nothing to do?
  if (!isFile(this)) return(naValue);

  # TO DO: Support gzipped files. /HB 2013-06-20
  if (isGzipped(this)) {
    return(naValue);
  }
  pathname <- getPathname(this);
  geometry <- memoizedCall2(this, function(this, ...) Biostrings::fastq.geometry(pathname));
  geometry;
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


setMethodS3("writeSample", "FastqDataFile", function(this, pathname, n, ordered=FALSE, ..., full=FALSE) {
  require("ShortRead") || throw("Package not loaded: ShortRead");

  # Argument 'pathname':
  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=TRUE);

  # Argument 'n':
  n <- Arguments$getInteger(n, range=c(1,Inf));

  # Argument 'ordered':
  ordered <- Arguments$getLogical(ordered);

  # Argument 'full':
  full <- Arguments$getLogical(full);


  pathnameFQ <- getPathname(this);

  # TODO: Added ram/buffer size option. /HB 2013-07-01
  fs <- FastqSampler(pathnameFQ, n=n, ordered=ordered, ...);
  on.exit({
    if (!is.null(fs)) close(fs);
  });
  data <- yield(fs);
  writeFastq(data, file=pathname, mode="w", full=full);
  data <- NULL;
  close(fs);
  fs <- NULL;

  newInstance(this, pathname);
}, protected=TRUE)


setMethodS3("findMateFile", "FastqDataFile", function(this, mustExist=FALSE, ...) {
  path <- getPath(this);
  filename <- getFilename(this);

  # Recognized R1/R2 filename patterns
  formats <- c(
    "_(%d)()",
    "_(%d)(_[0-9]+)",
    "_(R%d)()",
    "_(R%d)(_[0-9]+)"
  );
  formats <- sprintf("^(.*)%s[.]((fq|fastq)(|[.]gz))$", formats);

  # For R1 and R2...
  pathnameM <- NULL;
  for (mm in 1:2) {
    patterns <- unname(sapply(formats, FUN=sprintf, mm));
    pos <- unlist(sapply(patterns, FUN=regexpr, filename), use.names=FALSE);
    pattern <- patterns[pos != -1L][1L];
    if (!is.na(pattern)) {
      p1 <- gsub(pattern, "\\1", filename);
      p2 <- gsub(pattern, "\\2", filename);
      p3 <- gsub(pattern, "\\3", filename);
      ext <- gsub(pattern, "\\4", filename);
      mate <- 2L-mm+1L;
      p2M <-gsub(mm, mate, p2, fixed=TRUE);
      filenameM <- sprintf("%s_%s%s.%s", p1, p2M, p3, ext);
      pathnameM <- Arguments$getReadablePathname(filenameM, path=path, mustExist=FALSE);

      # Found mate file?
      if (isFile(pathnameM)) {
        break;
      }
    }
  } # for (mm ...)

  # Sanity check
  if (mustExist && is.null(pathnameM)) {
    throw("Failed to locate mate-pair file: ", getPathname(this));
  }

  pathnameM;
}, protected=TRUE)

setMethodS3("getMateFile", "FastqDataFile", function(this, ...) {
  pathnameM <- findMateFile(this, ..., mustExist=TRUE);
  newInstance(this, pathnameM);
}, protected=TRUE)


############################################################################
# HISTORY:
# 2014-01-07
# o Now readGeometry() for FastqDataFile returns missing values for
#   non-existing files.
# 2013-11-02
# o Now FastqDataFile(Set) handles more paired filename formats.
# 2013-08-24
# o Added getMateFile() for FastqDataFile.
# 2013-07-10
# o BUG FIX: writeSample() would give Error in UseMethod("getPathname").
# 2013-07-01
# o Added writeSample() for FastqDataFile.
# 2013-06-25
# o Added getDefaultFullName() for FastqDataFile so <fullname>.fastq.gz
#   is properly handled.  Should ideally handled by R.filesets.
# 2013-06-20
# o Now readGeometry() and getGeometry() returns c(NA,NA) for
#   gzipped files. This is better than an error.
# 2012-06-28
# o Created.
############################################################################
