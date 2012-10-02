#!/usr/bin/env Rscript

############################################################################
#
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "annotationData/organisms/Human/";
filename <- "human_g1k_v37.fasta";
fa <- FastaReferenceFile(filename, path=path);
print(fa);

# Data set
dataSet <- "AlbertsonD_2012-SCC";
platform <- "Generic";
path <- file.path("fastqData", dataSet, platform);
ds <- IlluminaFastqDataSet$byPath(path);
print(ds);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build BWA index set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBwaIndexSet(fa, method="is", verbose=-10);
print(is);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Process at most two FASTQ files
ds <- extract(ds, 1:min(3, length(ds)));

# In addition to SAM read group data inferred from the Illumina FASTQ
# files, manual set the library information for the whole data set.
setSamReadGroup(ds, SamReadGroup(LB="MPS-034"));

# BWA with BWA 'aln' options '-n 2' and '-q 40'.
alg <- BwaAlignment(ds, indexSet=is, n=2, q=40);
print(alg);

bs <- process(alg, verbose=-20);
print(bs);

# Display an example BAM file
bf <- getFile(bs, 1);
print(bf);


############################################################################
# HISTORY:
# 2012-09-25
# o Created.
############################################################################
