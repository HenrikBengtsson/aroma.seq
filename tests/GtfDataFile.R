library("aroma.seq")
setOption(aromaSettings, "devel/parallel", "none")

# Setup (writable) local data directory structure
setupExampleData()

organism <- "SaccharomycesCerevisiae"

gtf <- GtfDataFile$byOrganism(organism)
print(gtf)

names <- getSeqNames(gtf)
str(names)
print(gtf)

