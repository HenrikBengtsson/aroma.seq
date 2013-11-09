#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");
library("qrqc");
library("R.rsp");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);


dataSet <- "LambdaVirusExample";
organism <- "LambdaPhage";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data set
path <- file.path("fastqData", dataSet, organism);
ds <- FastqDataSet$byPath(path);

# Keep only the small FASTQ files
ds <- extract(ds, "reads");
print(ds);

pdfs <- report(ds, verbose=verbose);
print(pdfs);


############################################################################
# HISTORY:
# 2012-12-06
# o Created.
############################################################################
