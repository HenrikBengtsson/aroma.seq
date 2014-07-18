setConstructorS3("AromaUnitNucleotideCountsFile", function(...) {
  this <- extend(AromaUnitTabularBinaryFile(...), "AromaUnitNucleotideCountsFile");

  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})


setMethodS3("getAromaUncFile", "AromaUgpFile", function(this, ...) {
  chipType <- getChipType(this);
  AromaUnitNucleotideCountsFile$byChipType(chipType, nbrOfUnits=nbrOfUnits(this), ...);
})


setMethodS3("getFilenameExtension", "AromaUnitNucleotideCountsFile", function(static, ...) {
  "unc";
}, static=TRUE, protected=TRUE);


setMethodS3("getExtensionPattern", "AromaUnitNucleotideCountsFile", function(static, ...) {
  "[.](unc)$";
}, static=TRUE, protected=TRUE)


setMethodS3("getDefaultColumnNames", "AromaUnitNucleotideCountsFile", function(this, ...) {
  c("A", "C", "G", "T", "N");
})

setMethodS3("getGcContent", "AromaUnitNucleotideCountsFile", function(this, ..., fields=c("G", "C"), na.rm=TRUE) {
  data <- readDataFrame(this, ...);

  # Drop unknown nucleotides?
  if (na.rm) {
    data <- data[setdiff(colnames(data), "N")];
  }

  # GC count
  y <- Reduce("+", data[fields])

  # Total count
  total <- Reduce("+", data);

  # GC fraction
  y <- y / total;

  y;
})


setMethodS3("allocateFromUgp", "AromaUnitNucleotideCountsFile", function(static, ugp=NULL, createdOn=Sys.time(), createdBy=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'ugp':
  ugp <- Arguments$getInstanceOf(ugp, "AromaUgpFile");

  nbrOfCols <- length(getDefaultColumnNames(static));
  unc <- AromaUnitNucleotideCountsFile$allocateFromUnitAnnotationDataFile(ugp, ..., types=rep("integer", times=nbrOfCols), sizes=rep(4L, times=nbrOfCols));

  ftr <- readFooter(unc);
  ftr$sources <- list(
    binning = list(
      fullname=getFullName(ugp),
      filename=getFilename(ugp),
      filesize=getFileSize(ugp),
      checksum=getChecksum(ugp)
    )
  );
  ftr$createdOn <- createdOn;
  ftr$createdBy <- createdBy;
  writeFooter(unc, ftr);

  unc;
}, static=TRUE)


setMethodS3("importFromBSgenome", "AromaUnitNucleotideCountsFile", function(this, db, ugp=NULL, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'db':
  db <- Arguments$getInstanceOf(db, "BSgenome");

  # Argument 'ugp':
  if (!is.null(ugp)) {
    ugp <- Arguments$getInstanceOf(ugp, "AromaUgpFile");
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Importing nucleotide counts from BSGenome");


  if (is.null(ugp)) {
    verbose && enter(verbose, "Inferring UGP file from file footer");
    ftr <- readFooter(this);
    binSource <- ftr$sources$binning;
    stopifnot(is.list(binSource));
    fullname <- binSource$fullname;
    stopifnot(!is.null(fullname));
    ugp <- AromaUgpFile$byChipType(fullname, nbrOfUnits=nbrOfUnits(this));
    verbose && exit(verbose);
  }
  verbose && print(verbose, ugp);


  verbose && enter(verbose, "Setting up chromosome sequence labels from UGP");
  chromosomes <- getChromosomes(ugp);
  labels <- sprintf("chr%d", chromosomes);
  if (length(chromosomes) >= 23) labels[23] <- "chrX";
  if (length(chromosomes) >= 24) labels[24] <- "chrY";
  if (length(chromosomes) >= 25) labels[25] <- "chrM";
  names(chromosomes) <- labels;
  verbose && print(verbose, chromosomes);

  verbose && enter(verbose, "Asserting that all chromosomes exists and sequences in the BSgenome source");
  # Sanity check
  stopifnot(all(is.element(names(chromosomes), seqnames(db))));
  verbose && exit(verbose);

  verbose && exit(verbose);


  verbose && enter(verbose, "Binning chromosome by chromosome");

  by <- NULL;
  for (kk in seq_along(chromosomes)) {
    chr <- chromosomes[kk];
    chrLabel <- names(chromosomes)[kk];
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", kk, chrLabel, length(chromosomes)));

    units <- getUnitsOnChromosome(ugp, chromosome=chr);
    verbose && cat(verbose, "Units on this chromosome:");
    verbose && str(verbose, units);

    # Nothing todo?
    if (length(units) == 0L) {
      verbose && cat(verbose, "No units. Skipping.");
      verbose && exit(verbose);
      next;
    }

    xOut <- getPositions(ugp, units=units);
    verbose && cat(verbose, "Number of bins: ", length(xOut));
    if (is.null(by)) {
      by <- median(diff(sort(xOut)), na.rm=TRUE);
      verbose && cat(verbose, "Inferred 'by': ", by);
    }
    bx <- c(xOut[1]-by/2, xOut+by/2);

    verbose && cat(verbose, "Nucleotide sequence:");
    seq <- db[[chrLabel]];
    verbose && print(verbose, seq);

    verbose && enter(verbose, "Binned counting of nucleotides");
    countsKK <- binTabulate(seq, bx=bx);
    verbose && exit(verbose);

    this[units,] <- countsKK;

    # Not needed anymore
    rm(seq, bx, xOut, countsKK, units);

    verbose && exit(verbose);
  } # for (kk ...)

  verbose && exit(verbose);


  verbose && enter(verbose, "Updating file footer");

  pkg <- db@seqs_pkgname;

  ftr <- readFooter(this);
  ftr$sources <- list(
    sequence = list(
      package=pkg,
      version=packageVersion(pkg)
    ),
    binning = list(
      fullname=getFullName(ugp),
      filename=getFilename(ugp),
      filesize=getFileSize(ugp),
      checksum=getChecksum(ugp)
    )
  );
  ftr$description <- sprintf("Nucleotides (A,C,G,T) according to the sequence source (%s) were counted for every bin defined by the UGP (%s). Positions with missing/unknown nucleotides where not counted.", pkg, getFilename(ugp));
  writeFooter(this, ftr);

  verbose && exit(verbose);


  verbose && exit(verbose);

  this;
}) # importFromBSgenome()



setConstructorS3("AromaUncFile", function(...) {
  this <- extend(AromaUnitNucleotideCountsFile(...), "AromaUncFile");
  # Parse attributes (all subclasses must call this in the constructor).
  if (!is.null(this$.pathname))
    setAttributesByTags(this);

  this;
})


############################################################################
# HISTORY:
# 2014-06-24
# o Added argument 'na.rm' to getGcContent().
# o Now this file type also allocates and 'N' column.
# 2012-10-18
# o Added getGcContent().
# 2012-10-16
# o Added getAromaUncFile() for AromaUgpFile.
# o Added allocateFromUgp().
# o Added importFromBSgenome().
# o Created from AromaUnitGcContentFile.R.
############################################################################
