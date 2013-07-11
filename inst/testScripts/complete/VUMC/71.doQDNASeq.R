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

# Process at most three FASTQ files
ds <- extract(ds, 1:min(3, length(ds)));

# In addition to SAM read group data inferred from the Illumina FASTQ
# files, manual set the library information for the whole data set.
setSamReadGroup(ds, SamReadGroup(LB="MPS-034", ID=1L));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Complete QDNASeq pipeline
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
cns <- doQDNASeq(ds, reference=fa, binWidth=1000, verbose=-20);
print(cns);


############################################################################
# HISTORY:
# 2013-07-11
# o Created from 21.BwaAlignment.R.
############################################################################
