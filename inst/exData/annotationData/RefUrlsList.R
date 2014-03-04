# source("RefUrlsList.R")

# Reference genomes and gene models from Ensembl for various organisms
RefUrlsList <- 
  list(
    
    # Cf. ftp://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/
    HS = # homo sapiens, ca. 201309
      c(reference="ftp://ftp.ensembl.org/pub/release-73/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.73.dna.primary_assembly.fa.gz",
        geneModel="ftp://ftp.ensembl.org/pub/release-73/gtf/homo_sapiens/Homo_sapiens.GRCh37.73.gtf.gz"),
    
    DM = # Drosophila Melanogaster, ca. 201305
      c(reference="ftp://ftp.ensembl.org/pub/release-70/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa.gz",
        geneModel="ftp://ftp.ensembl.org/pub/release-70/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.70.gtf.gz"),

    DM70 = # Used in Anders et al Nature Protocol
      c(reference="ftp://ftp.ensembl.org/pub/release-70/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa.gz",
        geneModel="ftp://ftp.ensembl.org/pub/release-70/gtf/drosophila_melanogaster/Drosophila_melanogaster.BDGP5.70.gtf.gz"),
    
    SC = # SaccharomycesCerevisiae, ca. 201307
      c(reference="ftp://ftp.ensembl.org/pub/release-73/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.EF4.73.dna.toplevel.fa.gz",
        geneModel="ftp://ftp.ensembl.org/pub/release-73/gtf/saccharomyces_cerevisiae/Saccharomyces_cerevisiae.EF4.73.gtf.gz")
  )
save(RefUrlsList, file="RefUrlsList_000.RData")



