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
# Build BWA index set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBwaIndexSet(fa, method="is", verbose=-10);
print(is);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BWA with BWA 'aln' options '-n 2' and '-q 40'.
alg <- BwaAlignment(ds, indexSet=is, n=2, q=40);
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
