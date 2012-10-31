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
dataSet <- "AlbertsonD_2012-Bladder_VUMC";
platform <- "Generic";
path <- file.path("bamData", dataSet, platform);
bs <- BamDataSet$byPath(path);
setFullNamesTranslator(bs, function(names, ...) {
  names <- gsub(".bowtie.sorted", "", names, fixed=TRUE);
  names <- gsub("([TN])$", ",\\1", names);
  names;
});
print(bs);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Count allele for chromosome 22
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chr <- 22;
units <- getUnitsOnChromosome(ugp, chr);
data <- readDataFrame(ugp, units=units);
data$chromosome <- sprintf("chr%d", data$chromosome);
filename <- sprintf("SNPs,chr%d.txt", chr);
res <- writeDataFrame(data, file=filename, col.names=FALSE, header=NULL);
print(res);


ds <- process(bc, verbose=verbose);
verbose && print(verbose, ds);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize GC content
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bgn <- BinnedGcNormalization(ds);
verbose && print(verbose, bgn);

dsN <- process(bgn, verbose=verbose);
verbose && print(verbose, dsN);

# DEBUG...
if (FALSE) {
  muN <- median(yN, na.rm=TRUE);
  rho <- muN / median(y, na.rm=TRUE);
  xlab <- "GC content fraction";
  ylab <- "Signals";
  gcLim <- c(0.25,0.75);
  Clim <- c(0,6);
  plot(gc, rho*y, pch=".", col="#999999", xlim=gcLim, ylim=c(0,6), xlab=xlab, ylab=ylab);
  points(gc, yN, pch=".");
  abline(h=muN, col="#cccccc");
  stext(side=3, pos=1, line=-1, cex=0.8, sprintf("n=%d", sum(is.finite(yN))));
}


############################################################################
# HISTORY:
# 2012-10-21
# o Now we can do CbsModel(dsT, ref="constant(2)").
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
