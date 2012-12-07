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
print(ds);

pdfs <- report(ds, verbose=verbose);
print(pdfs);


############################################################################
# HISTORY:
# 2012-12-06
# o Created.
############################################################################
