# Download iGenomes annotation data [1,2]:
# url="ftp://igenome:G3nom3s4u@ussd-ftp.illumina.com/Drosophila_melanogaster/Ensembl/BDGP5.25/Drosophila_melanogaster_Ensembl_BDGP5.25.tar.gz"
# wget $url
#
# Download FASTQ data set [3] [~1.6Gb]:
# url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE32nnn/GSE32038/suppl/GSE32038_simulated_fastq_files.tar.gz"
# wget $url
# gunzip GSE32038_simulated_fastq_files.tar.gz
# mkdir GSE32038/DROME/
# tar -xvf GSE32038_simulated_fastq_files.tar -C GSE32038/DROME/
#
# REFERENCES:
# [1] http://tophat.cbcb.umd.edu/igenomes.shtml
# [2] http://www.uniprot.org/taxonomy/7227
# [3] http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE32038

library("aroma.seq")
setOption(aromaSettings, "devel/parallel", "BiocParallel::BatchJobs")
library("BatchJobs")

organism <- "DROME"

path <- file.path("annotationData", "organisms", organism)
fa <- FastaReferenceFile("genome.fa", path=path)
print(fa)

gtf <- GenericDataFile("genes.gtf", path=file.path(path, "Genes"))
print(gtf)

path <- file.path("fastqData", "GSE32038", "DROME")
fqs <- FastqDataSet$byPath(path, paired=TRUE)
print(fqs)

# Assert that PE file pairs are recognized
pairs <- getFilePairs(fqs)
stopifnot(identical(dim(pairs), c(6L,2L)))

bams <- doTopHat2(fqs, reference=fa, transcripts=gtf, verbose=-100)
print(bams)

