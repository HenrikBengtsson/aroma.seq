############################################################################
# 
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Indexing a reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "annotationData/organisms/LambdaPhage";
pathname <- file.path(path, "lambda_virus.fa");
pathname <- Arguments$getReadablePathname(pathname);
res <- bwaIndex(pathname, "-a"="is");
print(res);

files <- list.files(path=path);
print(files);


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################
