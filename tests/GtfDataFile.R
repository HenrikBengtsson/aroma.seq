library("aroma.seq")
setOption("R.filesets/parallel", "none")

# Setup (writable) local data directory structure
setupExampleData()

organism <- "SaccharomycesCerevisiae"

gtf <- GtfDataFile$byOrganism(organism)
print(gtf)

names <- getSeqNames(gtf)
str(names)
print(gtf)

