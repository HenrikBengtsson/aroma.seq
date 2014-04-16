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
#   \item{struct}{A directory structure format.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Paired reads}{
#   There is a community/industry convention for paired-end runs:
#   "The reads are reported two FASTQ files, such that the n:th read in
#   the first file is mate-paired to the n:th read in the second file.
#   The read IDs must match." [1]
#
#   If argument \code{paired=TRUE}, the \code{files} arguments is assumed
#   to contain the "R1" files whereas the corresponding "R2" files are
#   implicit (inferred and located by matching the filenames).
# }
#
# \references{
#  [1] Simon Anders,
#      \emph{High-throughput sequencing: Alignment and related topic},
#      (38 slides), EMBL Heidelberg, June 2013.\cr
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("FastqDataSet", function(files=NULL, ..., paired=FALSE, struct=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'paired':
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  paired <- Arguments$getLogical(paired);

  this <- extend(GenericDataFileSet(files=files, ...), c("FastqDataSet", uses("AromaSeqDataFileSet")),
    .paired = paired
  );

  if (!is.null(struct)) directoryStructure(this) <- struct;

  this;
})


setMethodS3("as.character", "FastqDataSet", function(x, ...) {
  this <- x;
  s <- NextMethod("as.character");
  s <- c(s, sprintf("Is paired: %s", isPaired(this)));
  s;
}, protected=TRUE)


setMethodS3("byPath", "FastqDataSet", function(static, ..., recursive=FALSE, struct=NULL, paired=FALSE, pattern="[.](fq|fastq)(|[.]gz)$", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'paired':
  paired <- Arguments$getLogical(paired);

  # Argument 'recursive':
  recursive <- Arguments$getLogical(recursive);

  # Argument 'struct':
  if (!is.null(struct)) {
    struct <- directoryStructure(struct);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Setting up ", class(static)[1L], " by path");
  verbose && cat(verbose, "Recursive: ", recursive);
  verbose && cat(verbose, "Filename pattern: ", pattern);

  # SPEEDUP: This will make it only locate R1 FASTQ files
  if (paired && missing(pattern)) {
    patternT <- sprintf("(_1|_R1).*%s", pattern);
    verbose && cat(verbose, "Adjusted filename pattern to speedup setup for paired-end data set.");
  } else {
    patternT <- pattern;
  }
  verbose && cat(verbose, "Filename pattern: ", patternT);

  # Assume paired-end reads
  res <- NextMethod("byPath", recursive=recursive, pattern=patternT, paired=paired);

  # Paired reads?
  if (isPaired(res)) {
    verbose && enter(verbose, "Adjusting for paired-end reads");

    filenames <- basename(getPathnames(res));
    fullnames <- gsub(pattern, "", filenames);
    verbose && cat(verbose, "Full names: ", hpaste(fullnames));

    # Several alternatives exists:
    # (a) Does the fullnames end with _1 or _2?
    patterns <- c("_1$", "_R1$", "_R1(|_[0-9]+)$");
    for (pattern in patterns) {
      idxs <- grep(pattern, fullnames, fixed=FALSE);
      if (length(idxs) > 0L) break;
    }
    verbose && cat(verbose, "R1 filename pattern: ", pattern);

    if (length(idxs) == 0L) {
      throw("Failed to identify the R1 FASTQ files.");
    }

    res <- extract(res, idxs);

    verbose && exit(verbose);
  } # if (isPaired(res))

  verbose && exit(verbose);

  res;
}, protected=TRUE)


setMethodS3("getOrganism", "FastqDataSet", function(this, ...) {
  organism <- directoryItem(this, "organism");
  organism <- Arguments$getCharacter(organism, length=c(1L, 1L));
  organism;
}, protected=TRUE);


setMethodS3("isPaired", "FastqDataSet", function(this, ...) {
  this$.paired;
}, protected=TRUE);


setMethodS3("getFilePairs", "FastqDataSet", function(this, ...) {
  stopifnot(isPaired(this));

  pairs <- vector("list", length=2*length(this));
  dim(pairs) <- c(length(this), 2L);
  rownames(pairs) <- getNames(this);
  colnames(pairs) <- c("R1", "R2");

  for (ii in seq_along(this)) {
    r1 <- this[[ii]];
    if (isFile(r1)) {
      r2 <- getMateFile(r1);
    } else {
      r2 <- newInstance(r1, NA, mustExist=FALSE);
    }
    pairs[ii,1L] <- list(r1);
    pairs[ii,2L] <- list(r2);
  } # for (ii ...)

  pairs;
}, protected=TRUE)


setMethodS3("validate", "FastqDataSet", function(this, ...) {
  NextMethod("validate");
  if (isPaired(this)) {
    pairs <- getFilePairs(this);
  }
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


setMethodS3("findByName", "FastqDataSet", function(static, name, tags=NULL, organism=NULL, ..., paths="fastqData", pattern="[.](fq|fastq)(|[.]gz)$") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'organism':
  if (!is.null(organism)) {
    organism <- Arguments$getCharacter(organism);
  }

  NextMethod("findByPath", subdirs=organism, paths=paths, pattern=pattern);
}, static=TRUE)



setMethodS3("byName", "FastqDataSet", function(static, name, tags=NULL, organism=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'organism':
  if (!is.null(organism)) {
    organism <- Arguments$getCharacter(organism);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Setting up ", class(static)[1L], " by name");

  verbose && cat(verbose, "Organism: ", organism);

  suppressWarnings({
    paths <- findByName(static, name, tags=tags, organism=organism,
                                      firstOnly=FALSE, ...);
  })
  if (is.null(paths)) {
    path <- file.path(paste(c(name, tags), collapse=","), organism);
    throw("Cannot create ", class(static)[1], ".  No such directory: ", path);
  }

  verbose && cat(verbose, "Paths to possible data sets:");
  verbose && print(verbose, paths);

  # Record all exception
  exList <- list();

  res <- NULL;
  for (kk in seq_along(paths)) {
    path <- paths[kk];
    verbose && enter(verbose, sprintf("Trying path #%d of %d", kk, length(paths)));
    verbose && cat(verbose, "Path: ", path);

    tryCatch({
      suppressWarnings({
        res <- byPath(static, path=path, ..., verbose=verbose);
      });
    }, error = function(ex) {
      exList <<- append(exList, list(list(exception=ex, path=path)));

      verbose && cat(verbose, "Data set could not be setup for this path, because:");
      verbose && cat(verbose, ex$message);
    });

    if (!is.null(res)) {
      if (length(res) > 0) {
        verbose && cat(verbose, "Successful setup of data set.");
        verbose && exit(verbose);
        break;
      }
    }

    verbose && exit(verbose);
  } # for (kk ...)

  if (is.null(res)) {
    exMsgs <- sapply(exList, FUN=function(ex) {
      sprintf("%s (while trying '%s').",
                   ex$exception$message, ex$path);
    });
    exMsgs <- sprintf("(%d) %s", seq_along(exMsgs), exMsgs);
    exMsgs <- paste(exMsgs, collapse="  ");
    msg <- sprintf("Failed to setup a data set for any of %d data directories located. The following reasons were reported: %s", length(paths), exMsgs);
    verbose && cat(verbose, msg);
    throw(msg);
  }

  verbose && exit(verbose);

  res;
}, static=TRUE)





############################################################################
# HISTORY:
# 2014-01-07
# o Now getFilePairs() for FastqDataSet returns NA files (was R1) in
#   case mate file R2 could not be found.
# 2013-11-09
# o Added static findByName() and byName() for FastqDataSet.
# o Added getOrganism() to FastqDataSet.
# 2013-08-24
# o Added as.character() for FastqDataSet outputting also paired status.
# o Added byPath() for FastqDataSet that acknowledge 'paired' argument.
# o Added isPaired() and getFilePairs() for FastqDataSet.
# 2012-06-28
# o Created.
############################################################################
