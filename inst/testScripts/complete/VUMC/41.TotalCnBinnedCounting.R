#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");

dataSet <- "AlbertsonD_2012-Bladder_VUMC";
organism <- "HomoSapiens";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
by <- 50e3;
byTag <- sprintf("%dkb", by/1e3);
ugp <- AromaUgpFile$byChipType("GenericHuman", tags=byTag);
print(ugp);

# Data set
path <- file.path("bamData", dataSet, organism);
bs <- BamDataSet$byPath(path);
setFullNamesTranslator(bs, function(names, ...) {
  names <- gsub(".bowtie.sorted", "", names, fixed=TRUE);
  names <- gsub("([TN])$", ",\\1", names);
  names;
});
print(bs);


chrLabels <- names(getTargets(getFile(bs,1)));
chrMap <- c(1:25);
labels <- sprintf("chr%d", chrMap);
labels[23:25] <- sprintf("chr%s", c("X", "Y", "M"));
names(chrMap) <- labels;



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bin reads chromosome by chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bc <- TotalCnBinnedCounting(bs, targetUgp=ugp);
print(bc);

ds <- process(bc, verbose=verbose);
verbose && print(verbose, ds);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract tumor-normal pairs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsT <- extract(ds, sapply(bs, hasTag, "T"));
dsN <- extract(ds, sapply(bs, hasTag, "N"));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation and Chromosome Explorer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
seg <- CbsModel(dsT, ref=dsN, maxNAFraction=2/3);
verbose && print(verbose, seg);

ce <- ChromosomeExplorer(seg);
verbose && print(verbose, ce);
process(ce, verbose=verbose);



############################################################################
# HISTORY:
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
