############################################################################
# 
############################################################################
library("aroma.seq");


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reference genome
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
# Indexing of reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBwaIndexSet(fa, method="is", verbose=TRUE);
print(is);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reverse the order of the input set to test ordering
ds <- extract(ds, rev(seq(ds)));

alg <- BwaAlignment(ds, indexSet=is);
print(alg);

bs <- process(alg, verbose=TRUE);
print(bs);

# Sanity check
stopifnot(all(getFullNames(bs) == getFullNames(ds)));


############################################################################
# HISTORY:
# 2012-09-25
# o Can now run a complete BWA single-end alignment and get BAM files.
# 2012-09-24
# o Created.
############################################################################
