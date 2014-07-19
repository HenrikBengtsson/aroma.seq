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

# Order in which chromosomes are merged?
# (a) chr1, ..., chr9, chr10, ..., chr19, chrM, chrX, chrY
fas <- sortBy(fas, "mixedsort")
# (b) chr1, chr10, chr11, ..., chr19, chr2, ..., chr9, chrM, chrX, chrY
fas <- sortBy(fas, "lexicographic")
# (c) chr10, chr11, ..., chr19, chr1, chr2, ..., chr9, chrM, chrX, chrY
#     This is how it's done in ftp://ussd-ftp.illumina.com/Mus_musculus/
#     UCSC/mm9/Mus_musculus_UCSC_mm9.tar.gz
#     Doing it this way creates an identical file
fas <- sortBy(fas, "lexicographic")
fas <- extract(fas, order(nchar(getNames(fas)), decreasing=TRUE))

# Merge FASTA files
fa <- writeFastaReferenceFile(fas, filename="MusMusculus_mm9_chr1-22.fa", verbose=TRUE)
print(fa)


############################################################################
# HISTORY:
# 2014-02-20 [HB]
# o Created.
############################################################################
