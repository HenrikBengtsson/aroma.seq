############################################################################
#
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
unc <- AromaUncFile$allocateFromUgp(ugp, tags=c("HB20121016"), createdBy="Henrik Bengtsson, hb@biostat.ucsf.edu");
verbose && print(verbose, unc);

unc <- importFromBSgenome(unc, db=db, verbose=verbose);
verbose && print(verbose, unc);


# Calculate GC fractions for each bin
gc <- (Y[,"G"] + Y[,"C"]) / rowSums(Y, na.rm=TRUE);


############################################################################
# HISTORY:
# 2012-10-16
# o Created.
############################################################################
