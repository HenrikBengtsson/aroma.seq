############################################################################
#
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Indexing a reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "annotationData/organisms/Lambda_phage";
pathnameFA <- file.path(path, "lambda_virus.fa");
res <- bwaIndex(pathnameFA, method="bwtsw", verbose=TRUE);
print(res);

prefix <- bwaIndexPrefix(pathnameFA, method="bwtsw")
files <- list.files(path=getParent(prefix));
print(files);


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################
