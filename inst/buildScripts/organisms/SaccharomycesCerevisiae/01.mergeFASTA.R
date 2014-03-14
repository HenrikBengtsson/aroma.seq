# REFERENCES:
# ftp://hgdownload.cse.ucsc.edu/goldenPath/sacCer2/chromosomes/
library("aroma.seq")

organism <- "SaccharomycesCerevisiae"
ucscLabel <- "sacCer2"

# Setup all FASTA reference files
path <- dirname(FastaReferenceFile$findByOrganism(organism))
fas <- FastaReferenceSet$byPath(path)

# Reorder chrI,...,chrXVI,chrM,2micron
names <- getNames(fas)
names <- gsub("^chr", "", names)
chrs <- utils:::.roman2numeric(names)
o <- order(chrs, na.last=TRUE)
fas <- fas[o]
print(fas)

# Sanity check
stopifnot(length(fas) == 18L)


# Merge FASTA files
filename <- sprintf("%s_%s_chr1-%d.fa", organism, ucscLabel, length(fas))
fa <- writeFastaReferenceFile(fas, filename=filename, verbose=TRUE)
print(fa)


############################################################################
# HISTORY:
# 2014-03-12 [HB]
# o Created.
############################################################################
