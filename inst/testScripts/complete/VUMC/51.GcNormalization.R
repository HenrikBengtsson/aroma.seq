#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
by <- 50e3;
byTag <- sprintf("%dkb", by/1e3);
ugp <- AromaUgpFile$byChipType("GenericHuman", tags=byTag);
print(ugp);

unc <- getAromaUncFile(ugp);
print(unc);

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
# Normalize GC content
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get GC contents for all units
ugp <- getAromaUgpFile(ds);
unc <- getAromaUncFile(ugp);
counts <- as.matrix(unc[]);
gc <- rowSums(counts[,c("G","C")])/rowSums(counts);


for (ii in seq(ds)) {
  df <- getFile(ds, ii);
  name <- getFullName(df);
  verbose && enter(verbose, sprintf("Sample %d ('%s') of %d", ii, name, length(ds)));

  y <- extractMatrix(df, drop=TRUE, verbose=less(verbose, 10));
  verbose && cat(verbose, "Signals:");
  verbose && str(verbose, y);

  verbose && enter(verbose, "Normalizing signals (on the log scale) for GC content");
  ly <- log2(y);
  targetFcn <- function(...) 1;
  lyN <- normalizeGcContent(ly, gcContent=gc, targetFcn=targetFcn, .isLogged=TRUE, .returnFit=TRUE);
  yN <- 2^lyN;
  fit <- attr(lyN, "modelFit");
  verbose && cat(verbose, "Model fit:");
  verbose && str(verbose, fit);
  verbose && exit(verbose);


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

  verbose && exit(verbose);
} # for (ii ...)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Extract tumor-normal pairs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dsT <- extract(ds, sapply(bs, hasTag, "T"));
dsN <- extract(ds, sapply(bs, hasTag, "N"));

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation and Chromosome Explorer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
seg <- CbsModel(dsT, ref=dsN);
verbose && print(verbose, seg);

ce <- ChromosomeExplorer(seg);
verbose && print(verbose, ce);
process(ce, maxNAFraction=2/3, verbose=verbose);



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
