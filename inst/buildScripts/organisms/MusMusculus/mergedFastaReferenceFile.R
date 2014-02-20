# REFERENCES:
# ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/chromosomes/

# Setup all FASTA reference files
fas <- FastaReferenceSet$byPath("mm9")
fas <- sortBy(fas, "mixedsort")

# Extract the "base" ones
names <- gsub(".fa", "", getNames(fas), fixed=TRUE)
idxs <- grep("^chr([1-9]{1,2}|[XYM])$", names)
fas <- extract(fas, idxs)
print(fas)

print(getFullNames(fas))
