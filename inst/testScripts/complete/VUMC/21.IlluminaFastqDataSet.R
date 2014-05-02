#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Data set
dataSet <- "AlbertsonD_2012-SCC";
organism <- "HomoSapiens";
path <- file.path("fastqData", dataSet, organism);
ds <- IlluminaFastqDataSet$byPath(path);
print(ds);


df <- ds[[1]];
print(df);

rg <- getSamReadGroup(df);
print(rg);


############################################################################
# HISTORY:
# 2012-09-30
# o Created.
############################################################################
