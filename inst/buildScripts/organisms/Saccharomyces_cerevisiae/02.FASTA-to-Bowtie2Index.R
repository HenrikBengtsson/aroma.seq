library("aroma.seq")

organism <- "Saccharomyces_cerevisiae"

fa <- FastaReferenceFile$byOrganism(organism)
print(fa)
faMD5 <- getChecksumFile(fa)
print(faMD5)

is <- buildBowtie2IndexSet(fa, verbose=TRUE)
print(is)

