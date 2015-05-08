library("aroma.seq")
setOption("R.filesets/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")

# Setup (writable) local data directory structure
setupExampleData()


dataSet <- "TopHat-example"
organism <- "Lambda_phage"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=FALSE)
print(fqs)

## Split a FASTQ file into four parts
fq <- fqs[[1]]
print(fq)
fqsS <- splitUp(fq, size=1/4, path="tmp-splitUp", verbose=TRUE)
print(fqsS)

## Split a gzipped FASTQ file into four parts
fqZ <- newInstance(fq, gzip(fq, remove=FALSE))
print(fqZ)
fqsZS <- splitUp(fqZ, size=1/4, path="tmp-splitUp-gzip", gzip=FALSE, verbose=TRUE)
print(fqsZS)

## Assert output regardless of gzipped input or not
stopifnot(identical(getChecksum(fqsZS), getChecksum(fqsS)))


## Cleanup
removeDirectory("tmp-splitUp-gzip", recursive=TRUE)
file.remove(getPathname(fqZ))
removeDirectory("tmp-splitUp", recursive=TRUE)

