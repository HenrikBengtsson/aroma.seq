## ac <- GatkAlleleCounting(bs, targetUgp=ugp, fa=fa);
## print(ac);
## dsC <- process(ac, verbose=verbose);
## print(dsC);

setConstructorS3("GatkAlleleCounting", function(dataSet=NULL, targetUgp=NULL, fa=NULL, dropEmpty=TRUE, ..., .reqSetClass="BamDataSet") {
  require("aroma.cn") || throw("Package not loaded: aroma.cn");

  # Argument 'targetUgp':
  if (!is.null(targetUgp)) {
    targetUgp <- Arguments$getInstanceOf(targetUgp, "AromaUgpFile");
  }

  # Argument 'fa':
  if (!is.null(fa)) {
    fa <- Arguments$getInstanceOf(fa, "FastaReferenceFile");
  }

  # Argument 'dropEmpty':
  dropEmpty <- Arguments$getLogical(dropEmpty);


  extend(AromaTransform(dataSet=dataSet, ..., .reqSetClass=.reqSetClass), "GatkAlleleCounting",
    .targetUgp = targetUgp,
    .fa = fa,
    .dropEmpty = dropEmpty
  );
}) # GatkAlleleCounting()


setMethodS3("getParameters", "GatkAlleleCounting", function(this, ...) {
  params <- list(
    targetUgp = this$.targetUgp,
    fa = this$.fa,
    dropEmpty = this$.dropEmpty
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

setMethodS3("getExpectedOutputFullnames", "GatkAlleleCounting", function(this, ...) {
  names <- NextMethod("getExpectedOutputFullNames");
  names <- paste(names, "alleleCounts", sep=",");
  names;
}, protected=TRUE)


setMethodS3("getOutputDataSet0", "GatkAlleleCounting", function(this, ...) {
  NextMethod("getOutputDataSet0", className="TabularTextFileSet", pattern=".*,allel[e]*Counts[.]txt$");
})


setMethodS3("process", "GatkAlleleCounting", function(this, ..., overwrite=FALSE, verbose=FALSE) {


  verbose && enter(verbose, "Counting alleles for known SNPs");
  bs <- getInputDataSet(this);
  verbose && print(verbose, bs);

  params <- getParameters(this);
  targetUgp <- params$targetUgp;
  fa <- params$fa;
  dropEmpty <- params$dropEmpty;

  verbose && print(verbose, targetUgp);
  verbose && print(verbose, fa);

  # Get output path
  pathD <- getPath(this);
  verbose && cat(verbose, "Output path: ", pathD);

  bedf <- NULL;

  for (ii in seq(bs)) {
    bf <- getFile(bs, ii);
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", ii, getName(bf), length(bs)));

    filename <- sprintf("%s,alleleCounts.txt", getFullName(bf));
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
      verbose && cat(verbose, "GATK system result: ", res);

      # (d) Cleanup GATK output file and resave
      verbose && enter(verbose, "Parsing and cleaning up the GATK output file");
      db <- TabularTextFile(pathnameDT);
      verbose && print(verbose, db);
    
      colClassPatterns <- c("(Locus|_base_counts)"="character");
      data <- readDataFrame(db, colClassPatterns=colClassPatterns);
      verbose && str(verbose, data);
    
      # Parse (chromosome, position)
      chr <- gsub(":.*", "", data$Locus);
      pos <- gsub(".*:", "", data$Locus);
      pos <- as.integer(pos);
      stopifnot(length(pos) == length(chr));

      # Parse (A,C,G,T) counts
      # NOTE: Some allele have non-zero "N" counts, which happens
      # when a read is aligned to a SNP, but its nucleotide at the
      # SNP position could not be called.  Because of this, we need
      # to keep the "N" column as well.
      counts <- data[[ncol(data)]];
      rm(data);
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
      stopifnot(nrow(counts) == length(chr));

      # Summaries
      countsT <- counts[,c("A", "C", "G", "T"),drop=FALSE];

      # Per-nucleotide coverage
      alleleCoverage <- colSums(countsT);
      alleleCoverage <- as.integer(alleleCoverage);
      names(alleleCoverage) <- colnames(counts);

      # Total coverage
      totalCoverage <- sum(alleleCoverage);
      verbose && printf(verbose, "Total number of reads (with called alleles) covering a SNP: %d\n", totalCoverage);

      # Per-SNP coverage
      coverage <- rowSums(countsT, na.rm=TRUE);
      coverage <- as.integer(coverage);
      tblCoverage <- table(coverage);
      verbose && cat(verbose, "Distribution of SNP coverages:");
      verbose && print(verbose, tblCoverage);

      # Identify homozygous and heterozygous SNPs (with coverage >= 2)
      isHom <- rowAnys(counts == coverage);
      isHom[(coverage < 2L)] <- NA; # Unknown
      isHet <- !isHom;
      nbrOfHoms <- sum(isHom, na.rm=TRUE);
      nbrOfHets <- sum(isHet, na.rm=TRUE);
      verbose && printf(verbose, "Number of homozygous SNPs (with coverage >= 2): %d\n", nbrOfHoms);
      verbose && printf(verbose, "Number of heterozygous SNPs (with coverage >= 2): %d\n", nbrOfHets);

      tblHetCoverage <- table(isHet);
      verbose && cat(verbose, "Distribution of heterozygous SNP coverages:");
      verbose && print(verbose, tblHetCoverage);

      # Drop SNPs with no coverage?
      if (dropEmpty) {
        keep <- (coverage > 0L);
        if (sum(keep) < length(keep)) {
          chr <- chr[keep];
          pos <- pos[keep];
          counts <- counts[keep,,drop=FALSE];
        }
        rm(keep);
      } #if (dropEmpty)


      header <- list(
        description = "GATK DepthOfCoverage Results",
        totalCoverage = totalCoverage,
        alleleCoverage = sprintf("%s:%d", names(alleleCoverage), alleleCoverage),
        nbrOfHoms = nbrOfHoms,
        nbrOfHets = nbrOfHets,
        tblCoverage = sprintf("%s:%d", names(tblCoverage), tblCoverage),
        tblHetCoverage = sprintf("%s:%d", names(tblHetCoverage), tblHetCoverage),
        dropEmpty = dropEmpty
      );

      # Updated allele count data
      data <- data.frame(chromosome=chr, position=pos, counts);
      rm(ns, chr, pos, counts);

      # Store
      pathnameDTT <- sprintf("%s.tmp", pathnameDT); # AD HOC
      pathnameDTT <- Arguments$getWritablePathname(pathnameDTT, mustNotExist=TRUE);

      pkg <- aroma.seq;
      createdBy <- sprintf("%s v%s (%s)", getName(pkg), getVersion(pkg), getDate(pkg));
      writeDataFrame(data, file=pathnameDTT, header=header, createdBy=createdBy);
      rm(data, header);
      file.remove(pathnameDT); # AD HOC
      pathnameDT <- popTemporaryFile(pathnameDTT);

      verbose && exit(verbose);
      
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

  colClassPatterns <- c("*"="integer", "chromosome"="character");
  data <- readDataFrame(db, colClassPatterns=colClassPatterns);

  data;
}, protected=TRUE) # readGatkCountFile()


setMethodS3("getCombineBy", "GatkAlleleCounting", function(this, ...) {
  combineBy <- function(dfList) {
    # Drop loci with zero coverage
    countKeys <- c("A", "C", "G", "T", "N");
    dfList <- lapply(dfList, FUN=function(df) {
      dfT <- df[,countKeys];
      counts <- Reduce("+", dfT);
      df[counts > 0L,];
    });
  
    # Stack
    df <- Reduce(rbind, dfList);
    rm(dfList);
  
    dfK <- df[,c("chromosome", "position")];
    counts <- as.matrix(df[,countKeys]);
    keys <- Reduce(function(...) paste(..., sep="\t"), dfK);
    dups <- duplicated(keys);
    
    # Unique (chromosome,position)
    dfU <- dfK[!dups,];
    keysU <- keys[!dups];
    countsU <- counts[!dups,,drop=FALSE];
  
    # Remainders
    keys <- keys[dups];
    counts <- counts[dups,,drop=FALSE];
  
    while(length(keys) > 0L) {
      dups <- duplicated(keys);
  
      # Unique (chromosome,position)
      keysT <- keys[!dups];
      countsT <- counts[!dups,,drop=FALSE];
  
      # Add to total counts
      idxs <- match(keysT, table=keys);
      stopifnot(all(is.finite(idxs)));
      countsU[idxs,] <- countsU[idxs,] + countsT;
      
      # Remainders
      keys <- keys[dups];
      counts <- counts[dups,,drop=FALSE];
    } # while()
  
    data <- cbind(dfU, countsU);
    data;
  } # combineBy()

  combineBy;
}, protected=TRUE, static=TRUE)



############################################################################
# HISTORY:
# 2012-11-09
# o Added argument 'dropEmpty' to GatkAlleleCounting, which will drop
#   all SNPs with zero coverage from the output files.
# 2012-11-08
# o Implemented getOutputDataSet().
# 2012-11-02
# o Now process() also cleans up the GATK output file by dropping
#   redundant columns and parsing the allele counts into integers.
# 2012-10-31
# o Added GatkAlleleCounting.
# o Created.
############################################################################
