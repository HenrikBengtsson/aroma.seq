############################################################################
# 
############################################################################
library("aroma.seq");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Session information
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathH <- Sys.getenv("BT2_HOME");
printf("BT2_HOME: %s\n", pathH);

# Sanity check
pathH <- Arguments$getReadablePath(pathH);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create link to bowtie2's example/ directory
path <- file.path(pathH, "example");
createLink(target=path);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Indexing a reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- systemBowtie2Build("example/reference/lambda_virus.fa", "lambda_virus");
print(res);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res <- systemBowtie2("-x"="lambda_virus", "-U"="example/reads/reads_1.fq", "-S"="eg1.sam");
print(res);
 

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Results
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- "eg1.sam";
bfr <- readLines(pathname);
print(head(bfr));


############################################################################
# HISTORY:
# 2012-08-31
# o Created.
############################################################################
