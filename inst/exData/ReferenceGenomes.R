# source("ReferenceGenomes.R")

# Reference genomes and gene models from Ensembl for various organisms
RefUrlsList  <- 
  list(
    DrosophilaMelanogaster =  # ca. 201305
      c(reference="ftp://ftp.ensembl.org/pub/release-70/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa.gz",
        geneModel="ftp://ftp.ensembl.org/pub/release-70/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.70.gtf.gz"),
    
    SaccharomycesCerevisiae = # ca. 201307
      c(reference="ftp://ftp.ensembl.org/pub/release-72/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.EF4.72.dna.toplevel.fa.gz",
        geneModel="ftp://ftp.ensembl.org/pub/release-72/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.EF4.72.gtf.gz")
  )

# save(RefUrlsList, file="RefUrlsList_000.RData")



