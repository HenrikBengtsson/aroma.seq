############################################################################
# 
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Indexing a reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "annotationData/organisms/LambdaPhage";
pathnameFA <- file.path(path, "lambda_virus.fa");
res <- bwaIndex(pathnameFA, a="is", verbose=TRUE);
print(res);

files <- list.files(path=getParent(bwaIndexPrefix(pathnameFA)));
print(files);


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################
