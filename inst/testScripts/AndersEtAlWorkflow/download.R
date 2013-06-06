## source("step00.Downloads.R")
## - Download R packages, sample fastq files and metadata, reference genome seq

############################################
## Download requisite packages
bDownloadPkgs <- FALSE
if (bDownloadPkgs) {
    source( "http://www.bioconductor.org/biocLite.R" )
    biocLite("BiocUpgrade")
    biocLite( c("ShortRead","DESeq", "edgeR") )
}

############################################
## Download files from SRA

## Re. "SraRunInfo.csv" below:
## From the Anders et al preprint (http://arxiv.org/pdf/1302.3685v2.pdf):
## http://www.ncbi.nlm.nih.gov/sra?term=SRP001537 (the entire experiment corresponding to GEO accession GSE18508), users can download a table of the metadata into acomma-separated tabular le \SraRunInfo.csv". To do this, click on \Send to:" (top right corner), select \File", select format \RunInfo" and click on \Create File".

## Further reading:
## Perl script to download Run and Sample XML docs; contains hardcoded HMP project workarounds
## - http://www.hmpdacc.org/doc/get_SRA_run_and_sample_xml.pl


## Get metadata
sri = read.csv("SraRunInfo.csv", stringsAsFactors=FALSE)
keep = grep("CG8144|Untreated-",sri$LibraryName)
sri = sri[keep,]


## Download sra files
fs = basename(sri$download_path)   ## [... NOTE THE UNDERSCORE HERE IS PARTICULARLY DANGEROUS FOR PORTABLE CODE! ...]
if (bDoAll) {
    for(i in 1:nrow(sri))
        download.file(sri$download_path[i], fs[i])
}
stopifnot( all(file.exists(fs)) )

## Convert to fastq, splitting paired-end files
if (bDoAll) {
    for(f in fs) {
        cmd = paste("fastq-dump --split-3", f)
        cat(cmd,"\n")
        system(cmd)
    }
}

############################################
## Download reference genome
if (bDoAll) {
    cmds <- c("wget ftp://ftp.ensembl.org/pub/release-70/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa.gz",
              "gunzip Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa.gz",
              "wget ftp://ftp.ensembl.org/pub/release-70/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.70.gtf.gz",
              "gunzip Drosophila_melanogaster.BDGP5.70.gtf.gz")
    sapply(cmds, function(cmd)
       {
           cat(cmd, "\n")
           system(cmd)
       })
}

############################################
