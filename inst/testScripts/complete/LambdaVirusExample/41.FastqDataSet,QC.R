#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");
library("qrqc");
library("R.rsp");
verbose <- Arguments$getVerbose(-10, timestamp=TRUE);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data set
dataSet <- "LambdaVirusExample";
platform <- "Generic";
path <- file.path("fastqData", dataSet, platform);
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
