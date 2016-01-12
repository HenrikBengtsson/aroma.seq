library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
if (fullTest) {

message("*** BwaIndexSet ...")

## Setup (writable) local data directory structure
setupExampleData()

organisms <- c("Lambda_phage", "Saccharomyces_cerevisiae")

for (organism in organisms) {
  message("Organism: ", organism)
  fa <- FastaReferenceFile$byOrganism(organism)
  print(fa)

  message("BWA index set ...")
  ## Build index set
  is <- buildBwaIndexSet(fa, verbose=-10)
  print(is)
  stopifnot(isCompatibleWith(is, fa))
  stopifnot(isCompatibleWith(fa, is))
  message("BWA index set ... DONE")

  ## Build index set
  message("Bowtie2 index set ...")
  is <- buildBowtie2IndexSet(fa, verbose=-100)
  print(is)
  stopifnot(isCompatibleWith(is, fa))
  stopifnot(isCompatibleWith(fa, is))
  message("Bowtie2 index set ... DONE")
}

message("*** BwaIndexSet ... DONE")

} # if (fullTest)
