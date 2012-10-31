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
chr <- 22;
pathnameL <- sprintf("SNPs,chr%d.bed", chr);
if (!isFile(pathnameL)) {
  units <- getUnitsOnChromosome(ugp, chr);
  data <- readDataFrame(ugp, units=units);
  pathnameL <- writeDataFrame(data, file=pathnameL, col.names=FALSE, header=NULL);
  print(pathnameL);
}

bf <- getFile(bs, 1L);
print(bf);

pathnameD <- sprintf("%s,chr%d,allelCounts.txt", getFullName(bf), chr);
if (!isFile(pathnameD)) {
  res <- systemGATK(T="DepthOfCoverage", I=getPathname(bf), R=getPathname(fa), L=pathnameL,
         "--omitIntervalStatistics", "--omitLocusTable","--omitPerSampleStats", "--printBaseCounts",
         "o"=pathnameD, verbose=verbose);
}


############################################################################
# HISTORY:
# 2012-10-31
# o Created.
############################################################################
