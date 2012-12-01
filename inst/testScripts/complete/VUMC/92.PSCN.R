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
tags <- c();
#tags <- c(tags, "AB042");
platform <- "Generic";
path <- file.path("fastqData", dataSet, platform);
ds <- IlluminaFastqDataSet$byPath(path);
print(ds);

# Work with a subset of all FASTQ files?
if (!is.null(arrays)) {
  ds <- extract(ds, arrays);
  print(ds);
} else {
  arrays <- seq_along(ds);
}

arraysTag <- sprintf("arrays=%s", seqToHumanReadable(arrays));
print(arraysTag);


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
  colClassPatterns <- c("*"="integer", "chromosome"="character");
  cb <- GatkAlleleCounting$getCombineBy();
  dataAC <- readDataFrame(dsC, colClassPatterns=colClassPatterns, combineBy=cb, verbose=verbose);
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

fit(seg, verbose=verbose);

array <- 1L;
segList <- list();
chrs <- getChromosomes(seg);
chrLabels <- as.character(chrs);
chrLabels[chrs == 23] <- "X";
chrLabels[chrs == 24] <- "Y";
chrLabels[chrs == 25] <- "MT";
for (kk in seq_along(chrs)) {
  chr <- chrs[kk];
  chrLabel <- chrLabels[kk];
  chrTag <- sprintf("Chr%02d", chr);
  verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", kk, chrTag, length(chrs)));

  segsT <- getRegions(seg, arrays=array, chromosomes=chr, flat=TRUE, url=NULL, verbose=verbose);
  verbose && print(verbose, segsT);
  idxsT <- which(gp$chromosome == chrLabel);
  posT <- gp$position[idxsT];
  totalT <- total[idxsT];
  isHetT <- isHet[idxsT];
  
  verbose && printf(verbose, "Number of heterozygous SNPs: %d (out of %d)\n", sum(isHetT), length(idxsT));
  
  cat("Total SNP read counts:\n");
  verbose && print(verbose, table(totalT));
  
  # Only look at SNPs with 3 or more reads
  minTotal <- 3L;
  keepT <- (totalT >= minTotal);
  verbose && summary(verbose, keepT);
  
  verbose && printf(verbose, "Number of heterozygous SNPs with coverage >= %d: %d (out of %d)\n", minTotal, sum(isHetT[keepT]), length(idxsT[keepT]));
  
  # For each segment...
  mus <- sds <- double(nrow(segsT));
  snpCounts <- dhCounts <- integer(nrow(segsT));
  for (ss in seq_len(nrow(segsT))) {
    segsSS <- segsT[ss,];
    verbose && printf(verbose, "Segment #%d of %d\n", ss, nrow(segsT));
    verbose && print(verbose, segsSS);
    keepSS <- (segsSS$start <= posT & posT <= segsSS$stop);
    verbose && printf(verbose, "Number of SNPs in segment: %d\n", sum(keepSS));
    keepSS <- keepT & keepSS;
    verbose && printf(verbose, "Number of SNPs with coverage >= %d in segment: %d\n", minTotal, sum(keepSS));
    snpCounts[ss] <- sum(keepSS);
    dhCounts[ss] <- sum(keepSS & isHetT);

    idxsSS <- idxsT[keepSS];
    verbose && printf(verbose, "SNP indices:\n");
    verbose && str(verbose, idxsSS);

    rhoSS <- rho[idxsSS];
    mu <- mean(rhoSS, na.rm=TRUE);
    sd <- sd(rhoSS, na.rm=TRUE);
    verbose && printf(verbose, "Average DH in segment: %f\n", mu);
    mus[ss] <- mu;
    sds[ss] <- sd;
  } # for (ss ...)
  
  
  segsT$snpCounts <- snpCounts;
  segsT$dhCounts <- dhCounts;
  segsT$dhMean <- mus;
  segsT$dhSE <- sds / sqrt(dhCounts);
  verbose && print(verbose, segsT[,-1]);

  segList[[chrTag]] <- segsT;

  verbose && exit(verbose);
} # for (kk ...)

segsB <- Reduce(rbind, segList);
verbose && print(verbose, segsB[,-1]);

sampleName <- getNames(dsG)[1];
fullnameB <- paste(c(sampleName, arraysTag), collapse=",");
filenameB <- sprintf("%s.tsv.xls", fullnameB);
hdr <- list(
  srcFiles = getFullNames(dsG),
  arrays = arraysTag,
  minCoverage = minTotal
);
pathnameB <- writeDataFrame(segsB[,-1L], header=hdr, file=filenameB, overwrite=TRUE);


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
