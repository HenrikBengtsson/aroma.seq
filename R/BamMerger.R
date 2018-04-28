###########################################################################/**
# @RdocClass BamMerger
#
# @title "The BamMerger class"
#
# \description{
#  @classhierarchy
#
#  A BamMerger takes a @see "BamDataSet" as input and merges subsets of
#  @see "BamDataFile":s into single ones.  How the grouping is done, can
#  be controlled by a parameter.
# }
#
# @synopsis
#
# \arguments{
#  \item{dataSet}{An @see "BamDataSet".}
#  \item{groupBy}{A @character string or an explicit named @list,
#   specifying which input files should be processed together.}
#  \item{...}{Not used.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Supported operating systems}{
#   This method is available on Linux, OSX, and Windows [1].
# }
#
# @author HB
#
# \seealso{
#   Internally @see "Rsamtools::mergeBam" is used for merging, sorting,
#   and indexing.
# }
#*/###########################################################################
setConstructorS3("BamMerger", function(dataSet=NULL, groupBy=NULL, ...) {
  if (!is.null(dataSet)) {
    # Argument 'groupBy':
    if (is.character(groupBy)) {
      groupBy <- match.arg(groupBy, choices=c("name"));
    } else if (is.list(groupBy)) {
      # Validated below
    } else {
      throw("Invalid argument 'groupBy': ", mode(groupBy));
    }
  }

  this <- extend(SamTransform(dataSet, groupBy=groupBy, ...), c("BamMerger", uses("FileGroupsInterface")));

  # Argument 'groupBy':
  if (is.list(groupBy)) {
    validateGroups(this, groups=groupBy);
  }

  this;
})


setMethodS3("getAsteriskTags", "BamMerger", function(this, ...) {
  # The 'groupBy' tag
  params <- getParameters(this);
  groupBy <- params$groupBy;
  if (is.character(groupBy)) {
    groupBy <- sprintf("by%s", capitalize(groupBy));
  } else {
    # Generate a short checksum
    groupBy <- getChecksum(groupBy, algo="crc32");
  }

  c("merged", groupBy);
}, protected=TRUE)

setMethodS3("getOutputDataSet", "BamMerger", function(this, onMissing=c("drop", "NA", "error"), ...) {
  # Argument 'onMissing':
  onMissing <- match.arg(onMissing);


  ## Find all existing output data files
  path <- getPath(this);
  bams <- BamDataSet$byPath(path, ...);

  ## Order according to grouped input data set
  names <- getGroupNames(this);
  bams <- extract(bams, names, onMissing=onMissing);

  bams;
}) # getOutputDataSet() for BamMerger



setMethodS3("process", "BamMerger", function(this, ..., skip=TRUE, force=FALSE, verbose=FALSE) {
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


  verbose && enter(verbose, "Merging BAM files");
  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  groups <- getGroups(this);
  verbose && printf(verbose, "Merging into %d groups: %s\n", length(groups), hpaste(names(groups)));
  verbose && str(verbose, head(groups));

  if (force) {
    todo <- seq_along(groups);
  } else {
    todo <- findFilesTodo(this, verbose=less(verbose, 1));
    # Already done?
    if (length(todo) == 0L) {
      verbose && cat(verbose, "Already done. Skipping.");
      res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));
      verbose && exit(verbose);
      return(invisible(res));
    }
  }

  verbose && cat(verbose, "Number of files: ", length(ds));
  verbose && cat(verbose, "Number of groups: ", length(groups));

  path <- getPath(this);
  verbose && cat(verbose, "Output path: ", path);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply merger to each group of BAM files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- listenv()

  IDXS <- groups[todo]
  for (gg in seq_along(IDXS)) {
    idxs <- IDXS[[gg]]
    dfList <- as.list(ds[idxs])

    group <- names(IDXS)[gg]
    if (is.null(group)) group <- "<unknown>"

    verbose && enter(verbose, sprintf("Merging BAM set #%d ('%s') of %d", gg, group, length(IDXS)))
    verbose && cat(verbose, "Number of BAM files to be merged: ", length(dfList))

    pathnames <- unlist(sapply(dfList, FUN=getPathname), use.names=FALSE)
    name <- names(dfList)[1L]
    filenameD <- sprintf("%s.bam", name)
    filenameD <- Arguments$getFilename(filenameD)
    pathnameBAM <- file.path(path, filenameD)
    pathnameBAI <- sprintf("%s.bai", pathnameBAM)
    verbose && cat(verbose, "Output BAM file: ", pathnameBAM)
    verbose && cat(verbose, "Output BAM index file: ", pathnameBAI)

    # Already done?
    if (!force && isFile(pathnameBAM)) {
      verbose && cat(verbose, "Already processed. Skipping.")
      res[[gg]] <- pathnameBAM
      verbose && exit(verbose)
      next
    }

    res[[gg]] %<-% {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # BEGIN: ATOMIC PROCESSING
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      verbose && enter(verbose, "Atomic processing")

      pathnameBAMT <- sprintf("%s.tmp", pathnameBAM)
      verbose && cat(verbose, "Temporary output BAM file: ", pathnameBAMT)
      pathnameBAIT <- sprintf("%s.bai", pathnameBAMT)

      # Special case; nothing to merge, just copy
      if (length(pathnames) == 1L) {
        verbose && enter(verbose, "Copying single BAM file (nothing to merge)")
        pathname <- pathnames[1L]
        copyFile(pathname, pathnameBAMT, copy.mode=FALSE, overwrite=force)

        # Copy existing index file or create index from scratch?
        bam <- BamDataFile(pathname)
        bai <- getIndexFile(bam)
        if (!is.null(bai)) {
          verbose && enter(verbose, "Copying existing BAM index file")
          copyFile(getPathname(bai), pathnameBAIT, copy.mode=FALSE, overwrite=force)
          verbose && exit(verbose)
        } else {
          verbose && enter(verbose, "Indexing BAM file")
          pathnameBAIT <- indexBam(pathnameBAMT)
          verbose && cat(verbose, "Generated index file: ", pathnameBAIT)
          file.rename(pathnameBAIT, pathnameBAI)
          verbose && exit(verbose)
        }

        verbose && exit(verbose)
      } else {
        pathnameBAMT2 <- mergeBam(pathnames, destination=pathnameBAMT, indexDestination=TRUE, overwrite=force)
        verbose && cat(verbose, "Generated BAM file: ", pathnameBAMT)
        .stop_if_not(pathnameBAMT2 == pathnameBAMT)
        file.rename(pathnameBAIT, pathnameBAI)
      }

      file.rename(pathnameBAMT, pathnameBAM)
      verbose && exit(verbose)
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # END: ATOMIC PROCESSING
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


      verbose && cat(verbose, "Created BAM file: ", pathnameBAM)
      verbose && cat(verbose, "Created BAM index file: ", pathnameBAI)
      # Sanity checks
      pathnameBAM <- Arguments$getReadablePathname(pathnameBAM)
      pathnameBAI <- Arguments$getReadablePathname(pathnameBAI)

      pathnameBAM
    } ## %<-%

    verbose && exit(verbose)
  } ## for (gg ...)

  res <- resolve(res)

  res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 50))

  verbose && exit(verbose)

  invisible(res)
}) # process()



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TO DROP
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("validateGroups", "BamMerger", function(this, groups, ...) {
  # Input data set
  ds <- getInputDataSet(this);
  nbrOfFiles <- length(ds);

  # Sanity checks
  idxs <- unlist(groups, use.names=FALSE);
  idxs <- Arguments$getIndices(idxs, max=nbrOfFiles);
  if (length(idxs) < nbrOfFiles) {
    throw("One or more input BAM files is not part of any group.");
  } else if (length(idxs) > nbrOfFiles) {
    throw("One or more input BAM files is part of more than one group.");
  }

  if (is.null(names(groups))) {
    throw("The list of groups does not have names.");
  }

  invisible(groups);
}, protected=TRUE)
