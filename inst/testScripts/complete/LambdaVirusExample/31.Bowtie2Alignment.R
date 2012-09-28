#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "annotationData/organisms/LambdaPhage";
fa <- FastaReferenceFile("lambda_virus.fa", path=path);
print(fa);

# Data set
dataSet <- "LambdaVirusExample";
platform <- "Generic";
path <- file.path("fastqData", dataSet, platform);
ds <- FastqDataSet$byPath(path);
print(ds);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build BWA index set
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
