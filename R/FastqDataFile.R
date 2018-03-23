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
#   \item{paired}{If @TRUE, the data set contains paired-end reads,
#     otherwise not.}
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
setConstructorS3("FastqDataFile", function(..., paired=FALSE) {
  extend(GenericDataFile(...), c("FastqDataFile", uses("AromaSeqDataFile")),
    .paired = paired
  );
})

setMethodS3("as.character", "FastqDataFile", function(x, ...) {
  this <- x;
  s <- NextMethod("as.character");
  s <- c(s, sprintf("Is paired: %s", isPaired(this)));
  n <- nbrOfSeqs(this, fast=TRUE);
  s <- c(s, sprintf("Number of sequences: %s", n));
  s <- c(s, sprintf("Common width of sequences: %d", getCommonSeqWidth(this, fast=TRUE)));
  s;
}, protected=TRUE)


setMethodS3("isPaired", "FastqDataFile", function(this, ...) {
  this$.paired;
}, protected=TRUE)


setMethodS3("getDefaultFullName", "FastqDataFile", function(this, paired=isPaired(this), ...) {
  name <- NextMethod("getDefaultFullName");

  # Drop filename extension
  name <- gsub("[.](fastq|fq)$", "", name, ignore.case=TRUE);

  # Drop paired-end suffixes, e.g. '_1'?
  if (paired) {
    name <- gsub("_(1|R1|2|R2)$", "", name, ignore.case=TRUE);
  }

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


setMethodS3("getGeometry", "FastqDataFile", function(this, force=FALSE, fast=FALSE, ...) {
  geometry <- this$.geometry;

  # If not cached, return NAs immediately?
  if (is.null(geometry) && fast) geometry <- c(NA_integer_, NA_integer_)

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

  pathname <- getPathname(this);
  geometry <- memoizedCall2(this, function(this, ...) Biostrings::fastq.geometry(pathname));

  geometry;
}, private=TRUE)


setMethodS3("writeSample", "FastqDataFile", function(this, pathname, n, ordered=FALSE, seed=NULL, ..., full=FALSE, verbose=FALSE) {
  use("ShortRead")

  # Argument 'pathname':
  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=TRUE);

  # Argument 'n':
  n <- Arguments$getInteger(n, range=c(1,Inf));

  # Argument 'ordered':
  ordered <- Arguments$getLogical(ordered);

  # Argument 'seed':
  if (!is.null(seed)) seed <- Arguments$getInteger(seed);

  # Argument 'full':
  full <- Arguments$getLogical(full);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set the random seed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(seed)) {
    verbose && enter(verbose, "Setting (temporary) random seed");
    oldRandomSeed <- NULL;
    if (exists(".Random.seed", mode="integer")) {
      oldRandomSeed <- get(".Random.seed", mode="integer");
    }
    on.exit({
      if (!is.null(oldRandomSeed)) {
        .Random.seed <<- oldRandomSeed;
      }
    }, add=TRUE);
    verbose && cat(verbose, "The random seed will be reset to its original state afterward.");
    verbose && cat(verbose, "Seed: ", seed);
    set.seed(seed);
    verbose && exit(verbose);
  }


  # TODO: Added ram/buffer size option. /HB 2013-07-01
  pathnameFQ <- getPathname(this);
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
  newInstance(this, pathnameM, paired=TRUE)
}, protected=TRUE)


setMethodS3("splitUp", "FastqDataFile", function(this, size, path=getPath(this), ..., gzip=isGzipped(this), maxFileSize=4e9, verbose=FALSE) {
  # Argument 'size':
  size <- Arguments$getNumeric(size, range=c(0,Inf))

  # Argument 'path':
  path <- Arguments$getWritablePath(path)
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enterf(verbose, "Splitting %s into multiple files", class(this)[1])
  verbose && print(verbose, this)
  verbose && printf(verbose, "Argument 'size': %g\n", size)
  nseqs <- nbrOfSeqs(this)
  if (size < 1) {
    nseqsPerFile <- ceiling(size*nseqs)
  } else {
    nseqsPerFile <- Arguments$getInteger(size, range=c(1,Inf))
  }
  verbose && printf(verbose, "Number of sequences/reads per file: %d\n", nseqsPerFile)
  nlinesPerFile <- 4L*nseqsPerFile
  verbose && printf(verbose, "Number of lines per file: %d\n", nlinesPerFile)
  
  nbrOfFiles <- ceiling(nseqs / nseqsPerFile)
  verbose && printf(verbose, "Number of files: %d\n", nbrOfFiles)
  if (nbrOfFiles == 1) {
    msg <- sprintf("Split of FASTQ file would result in a single file. Skipping.")
    verbose && cat(verbose, msg)
    warning(msg)
    return(this)
  } else if (nbrOfFiles > 9999) {
    msg <- sprintf("Split of FASTQ file would result in > 9999 files. Adjust 'size' argument: %g", size)
    verbose && cat(verbose, msg)
    throw(msg)
  }
  nbrOfBytesPerFile <- ceiling(getFileSize(this) / nbrOfFiles)
  verbose && printf(verbose, "Average number of bytes per file: %d\n", nbrOfBytesPerFile)
  if (nbrOfBytesPerFile > maxFileSize) {
    throw("Average file size exceeds maxFileSize=%d bytes: %d", maxFileSize, nbrOfBytesPerFile)
  }

  pathnameFQ <- getPathname(this)

  if (isPaired(this)) {
    fullname <- getDefaultFullName(this, paired=FALSE)
    fmt <- gsub("_(1|2|R1|R2)$", "_part%04d_\\1.fq", fullname)
  } else {
    fmt <- sprintf("%s_part%%04d.fq", getFullName(this))
  }
  filenames <- sprintf(fmt, seq_len(nbrOfFiles))
  pathnames <- file.path(path, filenames)

  clazz <- Class$forName(class(this)[1])
  
  ## Already done?
  isFile <- file_test("-f", pathnames)
  if (all(isFile)) {
    verbose && cat(verbose, "Already processed. Skipping.")
    res <- lapply(pathnames, FUN=clazz)
    res <- FastqDataSet(res)
    verbose && print(verbose, res)
    verbose && exit(verbose)
    return(res)
  }
  
  if (any(isFile)) {
    existing <- pathnames[isFile]
    throw(sprintf("Cannot split FASTQ file (%s). Some of the output files already exists: [%d] %s", length(existing), hpaste(existing)))
  }

  ## Split atomically
  pathnamesT <- sapply(pathnames, FUN=pushTemporaryFile)

  ## Gzip output files?
  outfile <- if (gzip) gzfile else file

  infile <- if (isGzipped(pathnameFQ)) gzfile else file
  con <- infile(pathnameFQ, open="r")
  on.exit({
    if (!is.null(con)) close(con)
  })

  nlines <- 0L
  for (kk in seq_len(nbrOfFiles)) {
    pathnameT <- pathnamesT[kk]
    verbose && enterf(verbose, "File #%d ('%s') of %d", kk, pathnameT, nbrOfFiles)

    ## Read chunk
    bfr <- readLines(con, n=nlinesPerFile, warn=FALSE, ok=TRUE)
    nbfr <- length(bfr)
    if (nbfr %% 4 != 0) {
      throw(sprintf("Read %d lines (after having read %d lines previously), which is not a multiple of four: %s", nbfr, nlines, pathnameFQ))
    }

    ## Write chunk
    local({
      conT <- outfile(pathnameT, open="wb")
      on.exit(close(conT))
      writeLines(bfr, con=conT)
    })
    
    bfr <- NULL ## Not need anymore
    
    nlines <- nlines + nbfr
    verbose && printf(verbose, "Total number of lines processed: %d (%g%%) of %d\n", nlines, 100*nlines/(4L*nseqs), 4L*nseqs)
    verbose && exit(verbose)
  } # for (kk ...)

  verbose && printf(verbose, "Read/wrote in total %d lines and generated %d files\n", nlines, nbrOfFiles)

  close(con)
  con <- NULL

  ## Rename temporary files
  pathnames <- sapply(pathnamesT, FUN=popTemporaryFile)
  
  res <- lapply(pathnames, FUN=clazz)
  res <- FastqDataSet(res)
  verbose && print(verbose, res)

  verbose && exit(verbose)

  res
})
