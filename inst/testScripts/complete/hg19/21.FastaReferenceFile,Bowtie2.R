############################################################################
#
############################################################################
library("aroma.seq");

organism <- "HomoSapiens";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Indexing a reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("annotationData", "organisms", organism);
filename <- "human_g1k_v37.fasta";
fa <- FastaReferenceFile(filename, path=path);
print(fa);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build BWA index set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
is <- buildBowtie2IndexSet(fa, method="is", verbose=-10);
print(is);


############################################################################
# HISTORY:
# 2013-06-27
# o Created (from ditto for BWA).
############################################################################
