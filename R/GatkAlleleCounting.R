## ac <- GatkAlleleCounting(bs, fa=fa, ugp=ugp);
## print(ac);
## dsC <- process(ac, verbose=verbose);
## print(dsC);

setConstructorS3("GatkAlleleCounting", function(dataSet=NULL, fa=NULL, ugp=NULL, ...) {
  # Argument 'dataSet':
  if (!is.null(dataSet)) {
    dataSet <- Arguments$getInstanceOf(dataSet, "BamDataSet");
  }

  # Argument 'fa':
  if (!is.null(fa)) {
    fa <- Arguments$getInstanceOf(fa, "FastaReferenceFile");
  }

  # Argument 'ugp':
  if (!is.null(ugp)) {
    ugp <- Arguments$getInstanceOf(ugp, "AromaUgpFile");
  }

  extend(Object(), "GatkAlleleCounting",
    .dataSet = dataSet,
    .fa = fa,
    .ugp = ugp
  );
}) # GatkAlleleCounting()


setMethodS3("getFastaReferenceFile", "GatkAlleleCounting", function(this, ...) {
  this$.fa;
})

setMethodS3("getAromaUgpFile", "GatkAlleleCounting", function(this, ...) {
  this$.ugp;
})

setMethodS3("getInputDataSet", "GatkAlleleCounting", function(this, ...) {
  dataSet <- this$.dataSet;
  dataSet;
})

setMethodS3("getPath", "GatkAlleleCounting", function(this, ...) {
  bs <- getInputDataSet(this, ...);
  dataSet <- getFullName(bs);
  chipType <- getChipType(bs);
  rootPath <- "gatkData";
  path <- file.path(rootPath, dataSet, chipType);
  path <- Arguments$getWritablePath(path);
  path;
})


setMethodS3("process", "GatkAlleleCounting", function(this, ..., overwrite=FALSE, verbose=FALSE) {
  ugp <- getAromaUgpFile(this);
  fa <- getFastaReferenceFile(this);
  bs <- getInputDataSet(this);

  verbose && enter(verbose, "Counting alleles for known SNPs");
  verbose && print(verbose, bs);
  verbose && print(verbose, ugp);
  verbose && print(verbose, fa);

  dataSet <- getFullName(bs);
  chipType <- getChipType(bs);

  # Get output path
  pathD <- getPath(this);

  bedf <- NULL;

  for (ii in seq(bs)) {
    bf <- getFile(bs, ii);
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", ii, getName(bf), length(bs)));

    filename <- sprintf("%s,allelCounts.txt", getFullName(bf));
    pathnameD <- Arguments$getReadablePathname(filename, path=pathD, mustExist=FALSE);

    if (overwrite || !isFile(pathnameD)) {
      pathnameD <- Arguments$getWritablePathname(pathnameD, mustNotExist=!overwrite);

      # (a) Instead of having GATK build a missing FAI index file, build it here
      buildIndex(fa, verbose=verbose);
      
      # (b) Get the BED file representation of UGP file
      if (is.null(bedf)) {
        bedf <- writeBedDataFile(ugp, chrMap=c(X=23, Y=24, MT=25), verbose=verbose);
        verbose && print(verbose, bedf);
      }

      # (c) Call GATK
      verbose && enter(verbose, "Calling GATK DepthOfCoverage");
      pathnameDT <- pushTemporaryFile(pathnameD);
      verbose && cat(verbose, "Writing to temporary file: ", pathnameDT);
      res <- systemGATK(T="DepthOfCoverage", I=getPathname(bf), R=getPathname(fa), L=getPathname(bedf), "--omitIntervalStatistics", "--omitLocusTable", "--omitPerSampleStats", "--printBaseCounts", "o"=pathnameDT, verbose=verbose);
      verbose && cat(verbose, "System result: ", res);
      pathnameD <- popTemporaryFile(pathnameDT);
      verbose && exit(verbose);
    } # if (overwrite || !isFile(...))
  
    # Parse GATK output file
    db <- TabularTextFile(pathnameD);
    verbose && print(verbose, db);
  
    verbose && exit(verbose);
  } # for (ii ...)

  res <- getOutputDataSet(this);
  verbose && print(verbose, res);

  verbose && exit(verbose);

  res;
}) # process()



setMethodS3("readGatkCountFile", "GatkAlleleCounting", function(this, array, ..., verbose=FALSE) {
  ds <- getOutputDataSet(this);
  array <- Arguments$getIndex(array, max=length(ds));

  df <- getFile(ds, array);
  print(df);

  # Parse GATK results
  pathname <- getPathname(df);
  db <- TabularTextFile(pathname);
  print(db);

  colClassPatterns <- c("(Locus|_base_counts)"="character");
  data <- readDataFrame(db, colClassPatterns=colClassPatterns);

  # Parse (chromosome,position)
  chr <- gsub(":.*", "", data$Locus);
  pos <- gsub(".*:", "", data$Locus);
  pos <- as.integer(pos);

  # Parse (A,C,G,T) counts
  counts <- data[[2L]];
  counts <- gsub("[ACGTN]:", ":", counts);
  counts <- gsub(" ", "", counts, fixed=TRUE);
  counts <- gsub("^:", "", counts);
  counts <- strsplit(counts, split=":", fixed=TRUE);
  ns <- sapply(counts, FUN=length);
  stopifnot(all(ns == 5L));
  counts <- unlist(counts, use.names=FALSE);
  counts <- as.integer(counts);
  counts <- matrix(counts, ncol=5L, byrow=TRUE);
  colnames(counts) <- c("A", "C", "G", "T", "N");
  # Sanity check
  stopifnot(all(counts[,"N"] == 0L));
  counts <- counts[, c("A", "C", "G", "T"), drop=FALSE];

  # Summarize results
  coverage <- rowSums(counts, na.rm=TRUE);
  coverage <- as.integer(coverage);
  isHom <- rowAnys(counts == coverage);
  isHet <- !isHom;
  tt1 <- colSums(counts);
  tt2 <- table(coverage);
  printf("Total number of reads covering a SNP: %d\n", sum(tt1));
  cat("Distribution of SNP coverages:");
  print(tt2);
  cat("Distribution of SNP alleles:");
  print(tt1);
  tt1b <- colSums(counts[isHet,,drop=FALSE]);
  tt2b <- table(coverage[isHet]);
  printf("Number of heterozygous SNPs: %d\n", sum(tt1));
  cat("Distribution of heterozygous SNP coverages:");
  print(tt2b);
  cat("Distribution of heterzygous SNP alleles:");
  print(tt1b);

  # Store counts
  data <- data.frame(chromosome=chr, position=pos, coverage=coverage);
  data <- cbind(data, counts);

  data;
}, protected=TRUE) # readGatkCountFile()


############################################################################
# HISTORY:
# 2012-10-31
# o Added GatkAlleleCounting.
# o Created.
############################################################################
