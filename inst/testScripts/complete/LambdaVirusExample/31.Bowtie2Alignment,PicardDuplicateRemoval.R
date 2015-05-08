#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");

dataSet <- "LambdaVirusExample";
organism <- "Lambda_phage";

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



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Remove duplicated reads using Picard
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dr <- PicardDuplicateRemoval(bs);
print(dr);

bsU <- process(dr, verbose=-20);
print(bsU);



############################################################################
# HISTORY:
# 2012-10-02
# o Created.
############################################################################
