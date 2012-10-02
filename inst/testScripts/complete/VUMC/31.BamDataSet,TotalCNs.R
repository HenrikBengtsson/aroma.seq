#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
# Extract tumor-normal pairs
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bsT <- extract(bs, sapply(bs, hasTag, "T"));
bsN <- extract(bs, sapply(bs, hasTag, "N"));

# Extract one pair
bfT <- getFile(bsT, 1);
bfN <- getFile(bsN, 1);
print(bfT);
print(bfN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Get read (start) position
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataT <- readReadPositions(bfT, verbose=-20);
str(dataT);

dataN <- readReadPositions(bfN, verbose=-20);
str(dataN);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bin reads for a chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chr <- "chr3";
xT <- subset(dataT, rname == chr, select="pos")$pos;
xN <- subset(dataN, rname == chr, select="pos")$pos;

xRange <- range(c(range(xT),range(xN)));
xOut <- seq(from=0, to=xRange[2], by=50e3);

count <- function(x, ...) length(x);
yTs <- binnedSmoothing(xT, x=xT, xOut=xOut, FUN=count, verbose=-10);
yNs <- binnedSmoothing(xN, x=xN, xOut=xOut, FUN=count, verbose=-10);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Calculate binned total copy numbers
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
C <- 2*yTs/yNs;
x <- xOut;

sampleName <- getName(bfT);
toPNG(sampleName, tags=c(chr, "50kb", "TCN"), width=1024, aspectRatio=0.3, {
  Clim <- c(0,5);
  plot(x/1e6,C, ylim=Clim, xlab="Position (Mb)", ylab="CN ratio");
  stext(side=3, pos=0, sampleName);
  stext(side=3, pos=1, chr);
  stext(side=4, pos=1, sprintf("(nT,nN)=(%d,%d)", length(xT), length(xN)));
});



############################################################################
# HISTORY:
# 2012-10-02
# o Created.
############################################################################
