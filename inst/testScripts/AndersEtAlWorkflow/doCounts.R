

load("samples.RData")
samples$countf <- paste(samples$LibraryName, "count", sep=".")
save(samples, file="samples.RData")  ## NB: Overwriting

gf <- GeneModelFile
samFiles <- paste0(samples$LibraryName, "_sn.sam")
outFiles <- samples$countf)
gfFile <- "Drosophila_melanogaster.BDGP5.70.gtf"
sapply(1:length(samFiles), function(i)
  {
    samFile <- samFiles[i]
    outFile <- outFiles[i]
    htseqCount(samFile, gfFile,
               optionsVec=c(s="no", a="10"),
               outFile)
    cat("Running htseq-count on sample ", i, "\n")  ## aroma verbose output can probably stand in for this
  })
cat("htseq-count done\n")


