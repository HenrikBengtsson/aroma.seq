# REFERENCES:
# ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/
library("aroma.seq")

organism <- "MusMusculus"

# Setup all FASTA reference files
path <- dirname(FastaReferenceFile$findByOrganism(organism))
fas <- FastaReferenceSet$byPath(path)

# Extract the "base" ones
fas <- extract(fas, "chr([0-9]{1,2}|[XYM])")
print(fas)
# Sanity check
stopifnot(length(fas) == 22L)

# Human-readable sorting?
fas <- sortBy(fas, "mixedsort")

# Merge FASTA files
fa <- writeFastaReferenceFile(fas, filename="mm9,chr1-22.fa", verbose=TRUE)
print(fa)


############################################################################
# HISTORY:
# 2014-02-20 [HB]
# o Created.
############################################################################
