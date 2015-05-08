############################################################################
#
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Indexing a reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "annotationData/organisms/Lambda_phage";
filename <- "lambda_virus.fa";
pathnameFA <- file.path(path,filename);
res <- bwaIndex(pathnameFA, method="bwtsw");
print(res);

prefix <- bwaIndexPrefix(pathnameFA, method="bwtsw");
print(prefix);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The FASTQ file to be aligned
dataSet <- "LambdaVirusExample";
organism <- "Lambda_phage";
path <- file.path("fastqData", dataSet, organism);
filename <- "reads_1.fq";
pathnameFQ <- Arguments$getReadablePathname(filename, path=path);

# BWA file
path <- file.path("bwaData", dataSet, organism);
filename <- "reads_1.sai";
pathnameSAI <- Arguments$getWritablePathname(filename, path=path);
res <- systemBWA("aln", f=pathnameSAI, prefix, pathnameFQ, verbose=TRUE);
print(res);

# SAM file
path <- file.path("bwaData", dataSet, organism);
filename <- "reads_1.sam";
pathnameSAM <- Arguments$getWritablePathname(filename, path=path);
res <- systemBWA("samse", f=pathnameSAM, prefix, pathnameSAI, pathnameFQ, verbose=TRUE);
print(res);

files <- list.files(path=path);
print(files);


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################
