############################################################################
# 
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Indexing a reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "annotationData/organisms/LambdaPhage";
filename <- "lambda_virus.fa";
pathnameFA <- Arguments$getReadablePathname(filename, path=path);
res <- bwaIndex(pathnameFA, "a"="is");
print(res);

prefix <- bwaIndexPrefix(pathnameFA);
print(prefix);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The FASTQ file to be aligned
dataSet <- "LambdaVirusExample";
platform <- "Generic";
path <- file.path("fastqData", dataSet, platform);
filename <- "reads_1.fq";
pathnameFQ <- Arguments$getReadablePathname(filename, path=path);

# BWA file
path <- file.path("bwaData", dataSet, platform);
filename <- "reads_1.sai";
pathnameSAI <- Arguments$getWritablePathname(filename, path=path);
res <- bwaAln(pathnameFQ, indexPrefix=prefix, pathnameD=pathnameSAI, verbose=TRUE);
print(res);

# SAM file
path <- file.path("bwaData", dataSet, platform);
filename <- "reads_1.sam";
pathnameSAM <- Arguments$getWritablePathname(filename, path=path);
res <- bwaSamse(pathnameSAI=pathnameSAI, pathnameFQ=pathnameFQ, indexPrefix=prefix, pathnameD=pathnameSAM, verbose=TRUE);
print(res);

files <- list.files(path=path);
print(files);

print(res);

files <- list.files(path=path);
print(files);


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################
