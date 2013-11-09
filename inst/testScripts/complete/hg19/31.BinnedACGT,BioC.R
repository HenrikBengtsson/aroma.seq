############################################################################
# REQUIREMENTS:
# annotationData/
#  chipTypes/
#   GenericHuman/
#    GenericHuman,50kb,HB20090503.ugp [1]
#
# REFERENCES:
# [1] http://aroma-project.org/data/annotationData/chipTypes/GenericHuman/
############################################################################
library("aroma.seq");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup the genome annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
by <- 50e3;
byTag <- sprintf("%dkb", by/1e3);
ugp <- AromaUgpFile$byChipType("GenericHuman", tags=byTag);
verbose && print(verbose, ugp);

pkg <- "BSgenome.Hsapiens.UCSC.hg19";
library(pkg, character.only=TRUE);
db <- get(pkg, mode="S4");
verbose && print(verbose, db);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Tabulate (A,C,G,T) in bins defined by the UGP
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
unc <- tryCatch({
  getAromaUncFile(ugp);
}, error=function(ex) { NULL });

# Build if missing
if (is.null(unc)) {
  unc <- AromaUncFile$allocateFromUgp(ugp, tags=c("HB20121021"), createdBy="Henrik Bengtsson, hb@biostat.ucsf.edu");
  verbose && print(verbose, unc);

  unc <- importFromBSgenome(unc, db=db, verbose=verbose);
}
verbose && print(verbose, unc);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Summaries
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
library("R.devices");
library("aroma.light");

toPNG(getFullName(unc), tags="GcContent", width=840, aspectRatio=0.5, {
  par(mar=c(4,4,1,1)+0.1);
  Y <- as.matrix(readDataFrame(unc));
  # Calculate GC fractions for each bin (=unit)
  gc <- (Y[,"G"] + Y[,"C"]) / rowSums(Y, na.rm=TRUE);
  plotDensity(gc, lwd=2, xlab="Fraction of GC nucleotides per bin");
  stext(side=3, pos=0, getFullName(unc));
  stext(side=3, pos=1, line=-1, cex=0.8, sprintf("n=%d", sum(is.finite(gc))));
})


############################################################################
# HISTORY:
# 2012-10-16
# o Created.
############################################################################
