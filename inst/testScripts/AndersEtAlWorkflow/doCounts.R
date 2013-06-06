## source("doCounts.R")

load("samples.RData")
samples$countf <- paste(samples$LibraryName, "count", sep=".")
save(samples, file="samples.RData")  ## NB: Overwriting

gf <- GeneModelFile

## [... RUN AROMA.SEQ WRAPPER 'HTSEQCOUNT' SAY, HERE TO REPLACE THE FOLLOWING ...]
if (FALSE) {
  cmds <- paste0("htseq-count -s no -a 10 ", samples$LibraryName, "_sn.sam ",
                 gf," > ", samples$countf)
  sapply(cmds, function(cmd)
         {
           cat("Running cmd:", cmd, "\n")
           system(cmd)
         }, simplify=FALSE)
  cat("done\n")
}




