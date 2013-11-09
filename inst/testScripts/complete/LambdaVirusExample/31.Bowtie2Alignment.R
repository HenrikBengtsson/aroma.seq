#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");

dataSet <- "LambdaVirusExample";
organism <- "LambdaPhage";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("annotationData", "organisms", organism);
fa <- FastaReferenceFile("lambda_virus.fa", path=path);
print(fa);

# Data set
path <- file.path("fastqData", dataSet, organism);
ds <- FastqDataSet$byPath(path);
print(ds);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build Bowtie2 index set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBowtie2IndexSet(fa, verbose=-10);
print(is);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Bowtie2 with bowtie2 ...
alg <- Bowtie2Alignment(ds, indexSet=is);
print(alg);

bs <- process(alg, verbose=-20);
print(bs);

# Display an example BAM file
bf <- getFile(bs, 1);
print(bf);


############################################################################
# HISTORY:
# 2012-09-27
# o Created.
############################################################################
