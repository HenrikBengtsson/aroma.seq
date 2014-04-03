## source("preProcess.R")
## - Should convert this to a function (will be roughly a stand-in for
##  level 3 code; i.e. following are hints into aroma-level design


source("download.R")


######################################################################
## Build reference genome index using bowtie2

library(aroma.seq)
setwd(DirRefFa)

## ASK HENRIK FOR THE PREFERRED PATTERN HERE
BowtieRefIndexPrefix <- file.path(Organism, IndexFilePrefix)
if (!file.exists(paste(BowtieRefIndexPrefix, ".1.bt2", sep="")))
  {
    bowtie2Build(refReads=RefFaFiles,
                 bowtieRefIndexPrefix=file.path("DM", "BDGP5"))
  }


############################################
## Quality assess fastq files

setwd(DirData)
if (bQAFastq) {
  library("ShortRead")
  fqQC = qa(dirPath=".", pattern=".fastq$", type="fastq")
  report(fqQC, type="html", dest="fastqQAreport")
}


############################################
## Gather sample metadata (experimental design)
## [... THIS IS A MANUAL STEP ...]
## Example below

## Here, sri is the SRA Run Info file
setwd(DirData)
sri$LibraryName = gsub("S2_DRSC_","",sri$LibraryName) # trim label
samples = unique(sri[,c("LibraryName","LibraryLayout")])
for(i in seq_len(nrow(samples))) {
    rw = (sri$LibraryName==samples$LibraryName[i])
    if(samples$LibraryLayout[i]=="PAIRED") {
        samples$fastq1[i] = paste0(sri$Run[rw],"_1.fastq",collapse=",")
        samples$fastq2[i] = paste0(sri$Run[rw],"_2.fastq",collapse=",")
    } else {
        samples$fastq1[i] = paste0(sri$Run[rw],".fastq",collapse=",")
        samples$fastq2[i] = ""
    }
}
samples$condition <- "CTL"
samples$condition[grep("RNAi",samples$LibraryName)] <- "KD"
samples$shortname <- paste( substr(samples$condition,1,2),
                           substr(samples$LibraryLayout,1,2),
                           seq_len(nrow(samples)), sep=".")




