############################################################################
# 
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Indexing a reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "annotationData/organisms/LambdaPhage";
filename <- "lambda_virus.fa";
faPathname <- file.path(path, filename);
res <- bwaIndex(filename, path=path, "-a"="is");
print(res);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# The FASTQ file to be aligned
dataSet <- "LambdaVirusExample";
platform <- "Generic";
path <- file.path("fastqData", dataSet, platform);
filename <- "reads_1.fq";
pathname <- file.path(path, filename);
pathname <- Arguments$getReadablePathname(filename, path=path);

# The aligne BWA file
pathD <- file.path("bwaData", dataSet, platform);
filenameD <- "reads_1.sai";
pathnameD <- Arguments$getWritablePathname(filenameD, path=pathD);

res <- systemBWA("aln", faPathname, pathname, stdout=pathnameD, stderr="foo.log");
print(res);

files <- list.files(path=pathD);
print(files);


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################
