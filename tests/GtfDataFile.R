library("aroma.seq")

# Setup (writable) local data directory structure
setupExampleData()

organism <- "Saccharomyces_cerevisiae"

gtf <- GtfDataFile$byOrganism(organism)
print(gtf)

names <- getSeqNames(gtf)
str(names)
print(gtf)


