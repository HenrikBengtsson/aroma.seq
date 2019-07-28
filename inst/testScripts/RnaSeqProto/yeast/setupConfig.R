# source("setupConfig.R")
# Set up run config:  paths, name, organism

setupConfig <- function(pathRef, pathData, datasetName, organism, bPairedEnd, bOverwrite=TRUE)
{
  
  config <- list()
  
  # Path to reference fasta and gtf files; these will be copied locally
  config$pathRef <- pathRef
  
  # Path to input fastq files; these will be copied (or symlinked) locally
  config$pathData <- pathData
  
  # Path to aroma.seq package data
  config$pathExData <- system.file("exData", package="aroma.seq")
  
  # Dataset metadata
  config$datasetName <- datasetName
  config$organism <- organism
  # config$bPairedEnd <- FALSE  # [ Assign PE status later, at alignment stage ]
  # config$qualityEncoding <- "illumina" # c("sanger", "solexa", "illumina"); cf. qrqc:ReadSeqFile
  
  config$bPairedEnd <- bPairedEnd
  config$bOverwrite <- bOverwrite
  
  return(config)
}
