############################################################################
#
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Session information
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Assert that BWA is available
pathname <- findBWA();
printf("BWA executable: %s\n", pathname);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Indexing a reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "annotationData/organisms/LambdaPhage";
pathnameFA <- file.path(path, "lambda_virus.fa");
pathnameFA <- Arguments$getReadablePathname(pathnameFA);
prefix <- bwaIndexPrefix(pathnameFA, method="is");
prefix <- Arguments$getWritablePathname(prefix) ## TAT addition; skip the mustNotExist=T
res <- systemBWA("index", a="is", p=prefix, pathnameFA, verbose=TRUE);
print(res);

files <- list.files(path=getParent(prefix));
print(files);


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################
