## source("doTopHat.R")

setwd(DirData)
load("samples.RData")
gf <- GeneModelFile
bowind <- file.path("DM", "BDGP5")
with(samples,
 {
     for (i in 1:nrow(samples))
     {
         optionsVec <- c(gf, LibraryName[i], 5)
         names(optionsVec) <- c("G", "o", "p")
         Reads1 <- unlist(strsplit(fastq1[i], ","))
         if (is.null(fastq2[i]) || nchar(fastq2[i])==0)
         {
             Reads2 <- NULL
         } else {
             Reads2 <- unlist(strsplit(fastq2[i], ","))
         }
         cat("Aligning sample", i, "\n")
         tophat(bowtieRefIndexPref=bowind, reads1=Reads1, reads2=Reads2, optionsVec=optionsVec)
     }
     cat("done\n")
 })

