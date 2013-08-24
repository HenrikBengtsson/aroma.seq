###########################################################################/**
# @RdocClass FastqDataSet
#
# @title "The FastqDataSet class"
#
# \description{
#  @classhierarchy
#
#  An FastqDataSet object represents a set of @see "FastqDataFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "FastqDataFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
#   \item{paired}{If @TRUE, the data set contains paired-end reads,
#     otherwise not.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Paired reads}{
#   If argument \code{paired=TRUE}, the \code{files} arguments is assumed
#   to contain the "R1" files whereas the corresponding "R2" files are
#   implicit (inferred and located by matching the filenames).
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("FastqDataSet", function(files=NULL, ..., paired=FALSE) {
  # Argument 'paired':
  paired <- Arguments$getLogical(paired);

  extend(GenericDataFileSet(files=files, ...), "FastqDataSet",
    .paired = paired
  );
})
setMethodS3("as.character", "FastqDataSet", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character");
  class <- class(s);

  s <- c(s, sprintf("Is paired: %s", isPaired(this)));

  class(s) <- class;
  s;
}, protected=TRUE)


setMethodS3("byPath", "FastqDataSet", function(static, ..., paired=FALSE, pattern="[.](fq|fastq)(|[.]gz)$") {
  res <- NextMethod("byPath", pattern=pattern, ..., paired=paired);

  # Paired reads?
  if (isPaired(res)) {
    filenames <- basename(getPathnames(res));
    # Keep only R1 files (and assume corresponding R2 files exist)
    patternR1 <- gsub("^_(R1|1)", "", pattern);
    patternR1 <- sprintf("_(R1|1)%s", patternR1);
    idxs <- grep(patternR1, filenames, fixed=FALSE);
    res <- extract(res, idxs);
  } # for (ii ...)

  res;
}, protected=TRUE)


setMethodS3("isPaired", "FastqDataSet", function(this, ...) {
  this$.paired;
}, protected=TRUE);

setMethodS3("getFilePairs", "FastqDataSet", function(this, ...) {
  stopifnot(isPaired(this));

  pairs <- vector("list", length=2*length(this));
  dim(pairs) <- c(length(this), 2L);
  rownames(pairs) <- getNames(this);
  colnames(pairs) <- c("R1", "R2");

  names <- getNames(this);
  for (ii in seq_along(this)) {
    r1 <- getFile(this, ii);
    path <- getPath(r1);
    filenameR1 <- getFilename(r1);
    if (regexpr("_R1.", filenameR1, fixed=TRUE) != -1L) {
      filenameR2 <- gsub("_R1.", "_R2.", filenameR1, fixed=TRUE);
      names[ii] <- gsub("_R1", "", names[ii], fixed=TRUE);
    } else if (regexpr("_1.", filenameR1, fixed=TRUE) != -1L) {
      filenameR2 <- gsub("_1.", "_2.", filenameR1, fixed=TRUE);
      names[ii] <- gsub("_1", "", names[ii], fixed=TRUE);
    } else {
      throw("Unrecognized R1/R2 filename pattern: ", filenameR1);
    }
    pathnameR2 <- file.path(path, filenameR2);
    r2 <- newInstance(r1, pathnameR2);
    pairs[ii,1L] <- list(r1);
    pairs[ii,2L] <- list(r2);
  } # for (ii ...)

  rownames(pairs) <- names;

  pairs;
}, protected=TRUE)

setMethodS3("validate", "FastqDataSet", function(this, ...) {
  NextMethod("validate");
}, protected=TRUE)


setMethodS3("getDepth", "FastqDataSet", function(this, ...) {
  1L;
}, protected=TRUE);

setMethodS3("getDefaultSamReadGroup", "FastqDataSet", function(this, ...) {
  SamReadGroup();
})

setMethodS3("setSamReadGroup", "FastqDataSet", function(this, rg, ...) {
  # Argument 'rg':
  if (!is.null(rg)) {
    rg <- Arguments$getInstanceOf(rg, "SamReadGroup");
  }
  this$.rg <- rg;
  invisible(this);
})

setMethodS3("getSamReadGroup", "FastqDataSet", function(this, ...) {
  rg <- this$.rg;
  if (is.null(rg)) {
    rg <- getDefaultSamReadGroup(this, ...);
  }
  rg;
})


############################################################################
# HISTORY:
# 2013-08-24
# o Added as.character() for FastqDataSet outputting also paired status.
# o Added byPath() for FastqDataSet that acknowledge 'paired' argument.
# o Added isPaired() and getFilePairs() for FastqDataSet.
# 2012-06-28
# o Created.
############################################################################
