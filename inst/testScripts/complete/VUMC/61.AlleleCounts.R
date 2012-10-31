#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Annotation data
path <- "annotationData/organisms/Human/";
filename <- "human_g1k_v37.fasta";
fa <- FastaReferenceFile(filename, path=path);
print(fa);

tags <- c("SNPs", "chr1-25");
ugp <- AromaUgpFile$byChipType("GenericHuman", tags=tags);
print(ugp);

# Data set
dataSet <- "AlbertsonD_2012-SCC,bwa,is";
platform <- "Generic";
path <- file.path("bwaData", dataSet, platform);
bs <- BamDataSet$byPath(path);
print(bs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Count allele for chromosome 22
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## ac <- AlleleCounting(bs, snpUgp=ugp);
## print(ac);
## dsC <- process(ac, verbose=verbose);
## print(dsC);


# Validate
overwrite <- FALSE;
ugp <- Arguments$getInstanceOf(ugp, "AromaUgpFile");
fa <- Arguments$getInstanceOf(fa, "FastaReferenceFile");

verbose && enter(verbose, "Counting alleles for known SNPs");
verbose && print(verbose, bs);
verbose && print(verbose, ugp);

dataSet <- getFullName(bs);
chipType <- getChipType(bs);

# Get output path
rootPath <- "gatkData"
pathD <- file.path(rootPath, dataSet, chipType);

#chromosomes <- getChromosomes(ugp);
#chrTags <- sprintf("chr%s", seqToHumanReadable(chromosomes));

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

verbose && exit(verbose);




      # Parse GATK results
      db <- TabularTextFile(pathnameD);
      print(db);
      data <- readDataFrame(db, colClassPatterns=c("(Locus|_base_counts)"="character"));

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
#      data <- data.frame(chromosome=chr, position=pos, coverage=coverage);
#      data <- cbind(data, counts);



############################################################################
# HISTORY:
# 2012-10-31
# o Created.
############################################################################
