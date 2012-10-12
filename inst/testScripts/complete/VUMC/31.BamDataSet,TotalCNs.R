#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");
library("matrixStats");
library("GenomicRanges"); # GRanges()
library("IRanges"); # IRanges()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
by <- 50e3;
byTag <- sprintf("%dkb", by/1e3);
ugp <- AromaUgpFile$byChipType("GenericHuman", tags=byTag);
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

dsN <- process(bc, verbose=verbose);
verbose && print(verbose, dsN);

stop();

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
# Bin reads chromosome by chromosome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
for (kk in seq(along=chrLabels)) {
  chrLabel <- chrLabels[kk];
  verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", kk, chrLabel, length(chrLabels)));

  chr <- chrMap[chrLabel];
  units <- getUnitsOnChromosome(ugp, chromosome=chr);
  xOut <- getPositions(ugp, units=units);
  nOut <- length(xOut);
  verbose && printf(verbose, "Target loci: [%d] %s\n", nOut, hpaste(xOut));

  # Nothing to do?
  if (nOut == 0) {
    verbose && cat(verbose, "No target loci on chromosome. Skipping.");
    verbose && exit(verbose);
    next;
  }

  # Read (start) positions on current chromosome
  gr <- GRanges(seqnames=Rle(chrLabel), ranges=IRanges(-500e6, +500e6));
  xT <- readReadPositions(bfT, which=gr, verbose=-20)$pos;
  verbose && str(verbose, xT);
  xN <- readReadPositions(bfN, which=gr, verbose=-20)$pos;
  verbose && str(verbose, xN);

  # Bin
  bx <- c(xOut[1]-by/2, xOut+by/2);
  yTs <- binCounts(xT, bx=bx);
  yNs <- binCounts(xN, bx=bx);

  # Calculate binned total copy numbers (assuming diploid genome)
  C <- 2*yTs/yNs;
  x <- xOut;

  sampleName <- getName(bfT);
  chrTag <- sprintf("chr%02d", chr);
  toPNG(sampleName, tags=c(chrTag, byTag, "TCN"), width=1024, aspectRatio=0.3, {
    par(mar=c(5,4,2,2)+0.1);
    Clim <- c(0,5);
    plot(x/1e6,C, ylim=Clim, xlab="Position (Mb)", ylab="CN ratio");
    stext(side=3, pos=0, sampleName);
    stext(side=3, pos=1, chrTag);
    stext(side=3, pos=1, line=-1, cex=0.8, sprintf("n=%d @ %s", sum(is.finite(C)), byTag));
    stext(side=4, pos=0, cex=0.8, sprintf("nT/nN=%d/%d=%.2f", length(xT), length(xN), length(xT)/length(xN)));
  });
 
  verbose && exit(verbose);
} # for (kk ...)



############################################################################
# HISTORY:
# 2012-10-10
# o Now plotting whole-genome TCN tracks, where data is loaded chromosome
#   by chromosome. Also utilizing generic UGP files.
# 2012-10-02
# o Created.
############################################################################
