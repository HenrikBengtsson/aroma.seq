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
pathname <- file.path(path, "lambda_virus.fa");
pathname <- Arguments$getReadablePathname(pathname);
res <- systemBWA("index", "-a"="is", pathname);
print(res);

files <- list.files(path=path);
print(files);


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################
