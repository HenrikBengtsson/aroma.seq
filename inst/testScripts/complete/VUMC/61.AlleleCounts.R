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
ac <- AlleleCounting(bs, snpUgp=ugp);
print(ac);

dsC <- process(ac, verbose=verbose);
print(dsC);


chromosomes <- getChromosomes(ugp);
chrTags <- sprintf("chr%s", seqToHumanReadable(chromosomes));
for (kk in seq(along=chromosomes)) {
  chr <- chromosomes[kk];
  pathnameL <- sprintf("SNPs,chr%d.bed", chr);
  if (!isFile(pathnameL)) {
    units <- getUnitsOnChromosome(ugp, chr);
    data <- readDataFrame(ugp, units=units);
    dups <- duplicated(data);
    if (any(dups)) {
      printf("Dropping duplicated genome positions: %d\n", sum(dups));
      data <- data[!dups,];
      units <- units[!dups];
    }
    pathnameL <- writeDataFrame(data, file=pathnameL, col.names=FALSE, header=NULL);
    print(pathnameL);
  }

  for (ii in seq(bs)) {
    bf <- getFile(bs, ii);
    print(bf);

    pathnameD <- sprintf("%s,chr%d,allelCounts.txt", getFullName(bf), chr);
    if (!isFile(pathnameD)) {
      # Instead of having GATK build a missing FAI index file, build it here
      buildIndex(fa, verbose=verbose);
    
      res <- systemGATK(T="DepthOfCoverage", I=getPathname(bf), R=getPathname(fa), L=pathnameL,
             "--omitIntervalStatistics", "--omitLocusTable", "--omitPerSampleStats", "--printBaseCounts",
             "o"=pathnameD, verbose=verbose);

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
    }
  } # for (ii ...)
} # for (kk ...)





############################################################################
# HISTORY:
# 2012-10-31
# o Created.
############################################################################
