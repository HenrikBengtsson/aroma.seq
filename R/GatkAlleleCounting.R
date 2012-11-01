## ac <- GatkAlleleCounting(bs, targetUgp=ugp, fa=fa);
## print(ac);
## dsC <- process(ac, verbose=verbose);
## print(dsC);

setConstructorS3("GatkAlleleCounting", function(dataSet=NULL, targetUgp=NULL, fa=NULL, ..., .reqSetClass="BamDataSet") {
  require("aroma.cn") || throw("Package not loaded: aroma.cn");

  # Argument 'targetUgp':
  if (!is.null(targetUgp)) {
    targetUgp <- Arguments$getInstanceOf(targetUgp, "AromaUgpFile");
  }

  # Argument 'fa':
  if (!is.null(fa)) {
    fa <- Arguments$getInstanceOf(fa, "FastaReferenceFile");
  }


  extend(AromaTransform(dataSet=dataSet, ..., .reqSetClass=.reqSetClass), "GatkAlleleCounting",
    .targetUgp = targetUgp,
    .fa = fa
  );
}) # GatkAlleleCounting()


setMethodS3("getParameters", "GatkAlleleCounting", function(this, ...) {
  params <- list(
    targetUgp = this$.targetUgp,
    fa = this$.fa
  );
  params;
}, protected=TRUE);


setMethodS3("getRootPath", "GatkAlleleCounting", function(this, ...) {
  "gatkData";
}, private=TRUE)


setMethodS3("getPath", "GatkAlleleCounting", function(this, ...) {
  path <- NextMethod("getPath", create=FALSE);
  path <- dirname(path);
  params <- getParameters(this);
  targetUgp <- params$targetUgp;
  chipType <- getChipType(targetUgp, fullname=FALSE);

  # The full path
  path <- filePath(path, chipType);
  path <- Arguments$getWritablePath(path); 

  # Verify that it is not the same as the input path
  inPath <- getPath(getInputDataSet(this));
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  path;
})


setMethodS3("process", "GatkAlleleCounting", function(this, ..., overwrite=FALSE, verbose=FALSE) {


  verbose && enter(verbose, "Counting alleles for known SNPs");
  bs <- getInputDataSet(this);
  verbose && print(verbose, bs);

  params <- getParameters(this);
  targetUgp <- params$targetUgp;
  fa <- params$fa;

  verbose && print(verbose, targetUgp);
  verbose && print(verbose, fa);

  # Get output path
  pathD <- getPath(this);
  verbose && cat(verbose, "Output path: ", pathD);

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
        bedf <- writeBedDataFile(targetUgp, chrMap=c(X=23, Y=24, MT=25), verbose=verbose);
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
