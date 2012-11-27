#!/usr/bin/env Rscript

############################################################################
# DESCRIPTION:
# This script shows how to (1) align single-end DNAseq reads in FASTQ files
# to the human genome (FASTA file), (2) count the aligned reads (BAM files)
# in uniformely distributed 50kb bins (UGP file), (3) normalize the counts
# for amount of GC-content in each bins (UNC file), and (4) finally 
# segment the normalized DNAseq total copy-number counts using 
# Circular Binary Segmentation (CBS) and (5) generate an interactive
# Chromosome Explorer report viewable in the browser.
#
# REQUIREMENTS:
# fastqData/
#  AlbertsonD_2012-SCC/
#   Generic/
#    <sample>_<barcode>_L[0-9]{3}_R[12]_[0-9]{3}.fastq [private data]
#
# annotationData/
#  chipTypes/
#   GenericHuman/
#    GenericHuman,50kb,HB20090503.ugp [1]
#    GenericHuman,50kb,HB20121021.unc [1]
#  organisms/
#   Human/
#    human_g1k_v37.fasta [2]
#
# REFERENCES:
# [1] http://aroma-project.org/data/annotationData/chipTypes/GenericHuman/
# [2] ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/
#                                                   human_g1k_v37.fasta.gz
############################################################################
library("aroma.seq");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);
args <- commandArgs(asValues=TRUE, excludeReserved=TRUE)[-1L];
verbose && cat(verbose, "Command line arguments:");
verbose && str(verbose, args);

# Argument 'arrays':
arrays <- args$arrays;
if (!is.null(arrays)) {
  arrays <- sprintf("c(%s)", arrays);
  expr <- parse(text=arrays);
  arrays <- eval(expr);
  arrays <- Arguments$getIndices(arrays);
}
print(arrays)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reference genome
path <- "annotationData/organisms/Human/";
filename <- "human_g1k_v37.fasta";
fa <- FastaReferenceFile(filename, path=path);
print(fa);

# Binning reference (by 50kb)
ugp <- AromaUgpFile$byChipType("GenericHuman", tags="50kb");
print(ugp);

# Binned nucleotide-content reference
# (created from BSgenome.Hsapiens.UCSC.hg19 Bioconductor reference)
unc <- getAromaUncFile(ugp);
print(unc);

# FASTQ data set
dataSet <- "AlbertsonD_2012-SCC";
platform <- "Generic";
path <- file.path("fastqData", dataSet, platform);
ds <- IlluminaFastqDataSet$byPath(path);
print(ds);


# Work with a subset of all FASTQ files
if (!is.null(arrays)) {
  ds <- extract(ds, arrays);
  print(ds);
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build BWA index set (iff missing)
is <- buildBwaIndexSet(fa, method="is", verbose=verbose);
print(is);

# In addition to SAM read group data inferred from the Illumina FASTQ
# files, manual set the library information for the whole data set.
setSamReadGroup(ds, SamReadGroup(LB="MPS-034"));

# BWA with BWA 'aln' options '-n 2' and '-q 40'.
alg <- BwaAlignment(ds, indexSet=is, n=2, q=40);
print(alg);

bs <- process(alg, verbose=verbose);
print(bs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Count reads per bin
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bc <- TotalCnBinnedCounting(bs, targetUgp=ugp);
print(bc);

ds <- process(bc, verbose=verbose);
verbose && print(verbose, ds);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize GC content (and rescale to median=2 assuming diploid)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bgn <- BinnedGcNormalization(ds);
verbose && print(verbose, bgn);

dsG <- process(bgn, verbose=verbose);
verbose && print(verbose, dsG);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Count allele
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
tags <- c("SNPs", "chr1-25");
ugpSNPs <- AromaUgpFile$byChipType("GenericHuman", tags=tags);
print(ugpSNPs);

ac <- GatkAlleleCounting(bs, targetUgp=ugpSNPs, fa=fa);
print(ac);
dsC <- process(ac, verbose=verbose);
print(dsC);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate TCN and DH
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract SNP allele counts
if (!exists("dataAC")) {
  cb <- GatkAlleleCounting$getCombineBy();
  dataAC <- readDataFrame(dsC, combineBy=cb, verbose=verbose);
  str(dataAC);
}

# Calculate DH
gp <- dataAC[,c("chromosome", "position")];
counts <- as.matrix(dataAC[,c("A", "C", "G", "T")]);

# Calculate TCN
total <- rowSums(counts);

# Naive calling of heterozygous SNPs
max <- rowMaxs(counts);
isHet <- (max != total);

# Calculate DH
rho <- max/total;
rho[!isHet] <- NA;



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Non-paired PSCBS where segmentation is based on TCN only
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
seg <- CbsModel(dsG, ref="constant(2)", maxNAFraction=2/3);
verbose && print(verbose, seg);

chr <- 4L;
segs <- getRegions(seg, arrays=1, chromosomes=chr, flat=TRUE, url=NULL, verbose=verbose);
print(segs);

segsT <- subset(segs, chromosome == chr);
idxsT <- which(gp$chromosome == chr);
posT <- gp$position[idxsT];
totalT <- total[idxsT];
isHetT <- isHet[idxsT];

printf("Number of heterozygous SNPs: %d (out of %d)\n", sum(isHetT), length(idxsT));

cat("Total SNP read counts:\n");
print(table(totalT));

# Only look at SNPs with 3 or more reads
minTotal <- 4L;
keepT <- (totalT >= minTotal);
print(summary(keepT));

printf("Number of heterozygous SNPs with coverage >= %d: %d (out of %d)\n", minTotal, sum(isHetT[keepT]), length(idxsT[keepT]));

# For each segment...
mus <- sds <- double(nrow(segsT));
counts <- integer(nrow(segsT));
for (ss in seq_len(nrow(segsT))) {
  segsSS <- segsT[ss,];
  printf("Segment #%d of %d\n", ss, nrow(segsT));
  print(segsSS);
  keepSS <- (segsSS$start <= posT & posT <= segsSS$stop);
  printf("Number of SNPs in segment: %d\n", sum(keepSS));
  keepSS <- keepT & keepSS;
  printf("Number of SNPs with coverage >= %d in segment: %d\n", minTotal, sum(keepSS));
  idxsSS <- idxsT[keepSS];
  printf("SNP indices:\n");
  str(idxsSS);
  rhoSS <- rho[idxsSS];
  mu <- mean(rhoSS, na.rm=TRUE);
  sd <- sd(rhoSS, na.rm=TRUE);
  printf("Average DH in segment: %f\n", mu);
  mus[ss] <- mu;
  sds[ss] <- sd;
  if (length(idxsSS) > 0L) counts[ss] <- length(idxsSS);
} # for (ss ...)


segsT$dhCounts <- counts;
segsT$dhMean <- mus;
segsT$dhSE <- sds / sqrt(counts);
print(segsT[,-1]);



############################################################################
# HISTORY:
# 2012-11-26
# o Now calculating average DH per TCN segment.
# 2012-10-21
# o Now a self-contained script from FASTQ to segmentation.
# o Now we can do CbsModel(..., ref="constant(2)").
# 2012-10-11
# o Now generating a Chromosome Explorer report.
# o Added TotalCnBinnedCounting() which calculates bin counts centered
#   at target loci specified by an UGP annotation file and outputs an
#   AromaUnitTotalCnBinarySet data set of DNAseq bin counts.
# 2012-10-10
# o Now plotting whole-genome TCN tracks, where data is loaded chromosome
#   by chromosome. Also utilizing generic UGP files.
# 2012-10-02
# o Created.
############################################################################
