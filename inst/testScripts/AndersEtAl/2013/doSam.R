## source("doSam.R")
## [... SOME OF THE FOLLOWING SHOULD PROBABLY BE OPTIONAL ...]


load("samples.RData")
for(i in seq_len(nrow(samples))) {
    lib = samples$LibraryName[i]
    ob = file.path(lib, "accepted_hits.bam")
    ## copy file
    system(paste0("cp ",ob," ",lib,".bam"))
    ## sort by name
    system(paste0("samtools sort -n ",lib,".bam ",lib,"_sn"))
    ## convert to SAM for htseq-count
    system(paste0("samtools view -o ",lib,"_sn.sam ",lib,"_sn.bam"))
    ## sort by position
    system(paste0("samtools sort ",ob," ",lib,"_s"))  ### [ TAT: I.e., can't sort by name and by position at the same time! ]
    ## for IGV
    system(paste0("samtools index ",lib,"_s.bam"))
}



