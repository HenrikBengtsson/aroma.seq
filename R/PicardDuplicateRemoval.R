###########################################################################/**
# @RdocClass PicardDuplicateRemoval
#
# @title "The PicardDuplicateRemoval class"
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
#  \item{dataSet}{An @see "BamDataSet".}
#  \item{tags}{Additional tags for the output data sets.}
#  \item{ASSUME_SORTED, VALIDATION_STRINGENCY, ...}{
#    Additional Picard MarkDuplicates arguments.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("PicardDuplicateRemoval", function(dataSet=NULL, tags="*", ASSUME_SORTED=TRUE, VALIDATION_STRINGENCY="SILENT", ...) {
  # Validate arguments
  if (!is.null(dataSet)) {
    # Argument 'dataSet':
    dataSet <- Arguments$getInstanceOf(dataSet, "BamDataSet");
  } # if (!is.null(dataSet))

  # Arguments '...':
  args <- list(
    ASSUME_SORTED=ASSUME_SORTED,
    VALIDATION_STRINGENCY=VALIDATION_STRINGENCY,
    ...
  );


  this <- extend(Object(), c("PicardDuplicateRemoval", uses("ParametersInterface")),
    .ds = dataSet,
    .tags = tags,
    .args = args
  );

  setTags(this, tags);

  this;
})


setMethodS3("as.character", "PicardDuplicateRemoval", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  ds <- getInputDataSet(this);
  s <- c(s, "Input data set:");
  s <- c(s, as.character(ds));

  # Additional arguments
  s <- c(s, sprintf("Picard arguments: %s", getParametersAsString(this)));

  class(s) <- "GenericSummary";
  s;
}, protected=TRUE)


setMethodS3("getInputDataSet", "PicardDuplicateRemoval", function(this, ...) {
  this$.ds;
})

setMethodS3("getOptionalArguments", "PicardDuplicateRemoval", function(this, ...) {
  this$.args;
}, protected=TRUE)


setMethodS3("getParameters", "PicardDuplicateRemoval", function(this, ...) {
  params <- NextMethod("getParameters");
  params <- c(params, getOptionalArguments(this));
  params;
}, protected=TRUE)


setMethodS3("getAsteriskTags", "PicardDuplicateRemoval", function(this, collapse=NULL, ...) {
  tags <- c("-dups");

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, protected=TRUE)


setMethodS3("getName", "PicardDuplicateRemoval", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
}, protected=TRUE)


setMethodS3("getFlavor", "PicardDuplicateRemoval", function(this, ...) {
  this$.flavor;
}, protected=TRUE)


setMethodS3("getTags", "PicardDuplicateRemoval", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input data set
  ds <- getInputDataSet(this);
  tags <- getTags(ds, collapse=collapse);

  # Get class-specific tags
  tags <- c(tags, this$.tags);

  # Update default tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0) {
    tags <- NULL;
  }

  tags;
}, protected=TRUE)


setMethodS3("setTags", "PicardDuplicateRemoval", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }

  this$.tags <- tags;
}, protected=TRUE)


setMethodS3("getFullName", "PicardDuplicateRemoval", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "PicardDuplicateRemoval", function(this, ...) {
  # Use same root path as input data set
  ds <- getInputDataSet(this);
  path <- getPath(ds);
  path <- getParent(path, depth=2L);
  # Sanity check
  stopifnot(regexpr("Data$", path) != -1);
  path;
}, protected=TRUE)


setMethodS3("getPath", "PicardDuplicateRemoval", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Organism
  ds <- getInputDataSet(this);
  organism <- getOrganism(ds);

  # The full path
  path <- filePath(rootPath, fullname, organism);

  if (create) {
    path <- Arguments$getWritablePath(path);
  } else {
    path <- Arguments$getReadablePath(path, mustExist=FALSE);
  }

  # Verify that it is not the same as the input path
  inPath <- getPath(ds);
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  path;
}, protected=TRUE)


setMethodS3("nbrOfFiles", "PicardDuplicateRemoval", function(this, ...) {
  ds <- getInputDataSet(this);
  length(ds);
})


setMethodS3("getOutputDataSet", "PicardDuplicateRemoval", function(this, ...) {
  ## Find all existing output data files
  path <- getPath(this);
  res <- BamDataSet$byPath(path, ...);

  ## Keep only those samples that exists in the input data set
  ds <- getInputDataSet(this);
  res <- extract(res, getFullNames(ds), onMissing="drop");

  ## TODO: Assert completeness
  res;
})



setMethodS3("process", "PicardDuplicateRemoval", function(this, ..., skip=TRUE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Picard removal of duplicated reads");

  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  nbrOfFiles <- nbrOfFiles(this);
  verbose && cat(verbose, "Number of files: ", nbrOfFiles);

  params <- getParameters(this);
  verbose && cat(verbose, "Additional Picard MarkDuplicates arguments:");
  verbose && str(verbose, params);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds, FUN=function(df, params, path, ...., skip=TRUE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'df':
    df <- Arguments$getInstanceOf(df, "BamDataFile");

    # Argument 'skip':
    skip <- Arguments$getLogical(skip);

    # Argument 'path':
    path <- Arguments$getWritablePath(path);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    verbose && enter(verbose, "Picard MarkDuplicates on one sample");

    pathname <- getPathname(df);
    verbose && cat(verbose, "Input BAM pathname: ", pathname);

    # Output BAM file
    fullname <- getFullName(df);
    filename <- sprintf("%s.bam", fullname);
    pathnameD <- Arguments$getWritablePathname(filename, path=path);
    verbose && cat(verbose, "Output BAM pathname: ", pathnameD);
    pathnameDI <- gsub("[.]bam$", ".bai", pathnameD);
    verbose && cat(verbose, "Output BAM index pathname: ", pathnameDI);

    # Nothing to do?
    done <- (skip && isFile(pathnameD) && isFile(pathnameDI));
    if (done) {
      verbose && cat(verbose, "Already processed. Skipping");
    } else {
      verbose && print(verbose, df);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (a) Filter via Picard
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!isFile(pathnameD)) {
        verbose && enter(verbose, "Calling Picard MarkDuplicates");

        pathnameM <- gsub("[.]bam$", ",picard,MarkDuplicates,metrics.log", pathnameD);
        args <- list(
          "MarkDuplicates",
          INPUT=pathname,
          OUTPUT=pathnameD,
          METRICS_FILE=pathnameM,
          REMOVE_DUPLICATES=TRUE,
          ASSUME_SORTED=TRUE,            # TODO: Assert! /HB 2012-10-02
          VALIDATION_STRINGENCY="SILENT"
        );
        verbose && cat(verbose, "Arguments:");
        verbose && str(verbose, df);

        # Assert no overwrite
        stopifnot(getAbsolutePath(pathnameD) != getAbsolutePath(pathname));
        stopifnot(getAbsolutePath(pathnameM) != getAbsolutePath(pathnameD));

        args$verbose <- less(verbose, 20);
        res <- do.call(systemPicard, args);
        verbose && cat(verbose, "System call results:");
        verbose && print(verbose, res);

        verbose && exit(verbose);
      } # if (!isFile(pathnameD))

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (b) Generic BAM index
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      bf <- BamDataFile(pathnameD);
      verbose && print(verbose, bf);

      if (!hasIndex(bf)) {
        verbose && enter(verbose, "Creating BAM index");
        bfi <- buildIndex(bf, verbose=less(verbose, 10));
        verbose && exit(verbose);
      }
    } # if (done)

    verbose && exit(verbose);

    invisible(list(pathnameD=pathnameD, pathnameDI=pathnameDI));
  }, params=params, path=getPath(this), skip=skip, verbose=verbose) # dsApply()

  res <- getOutputDataSet(this);

  verbose && exit(verbose);

  res;
})



############################################################################
# HISTORY:
# 2013-09-03
# o Now process() for PicardDuplicateRemoval utilizes dsApply().
# 2012-11-26
# o BUG FIX: getOutputDataSet() would return a data set with "missing"
#   files, if not complete.  Now it only returns the existing files.
# 2012-10-02
# o Created.
############################################################################
