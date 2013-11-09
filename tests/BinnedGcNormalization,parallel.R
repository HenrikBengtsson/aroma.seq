############################################################################
# DESCRIPTION:
# This script shows how to (1) align single-end DNAseq reads in FASTQ files
# to the human genome (FASTA file), (2) count the aligned reads (BAM files)
# in uniformely distributed 50kb bins (UGP file), (3) normalize the counts
# for amount of GC-content in each bins (UNC file), and (4) finally
# segment the normalized DNAseq total copy-number counts using
# Circular Binary Segmentation (CBS) and (5) generate an interactive
# Chromosome Explorer report viewable in the browser.
#
# REQUIREMENTS:
# fastqData/
#  AlbertsonD_2012-SCC/
#   HomoSapiens/
#    <sample>_<barcode>_L[0-9]{3}_R[12]_[0-9]{3}.fastq [private data]
#
# annotationData/
#  chipTypes/
#   GenericHuman/
#    GenericHuman,50kb,HB20090503.ugp [1]
#    GenericHuman,50kb,HB20121021.unc [1]
#  organisms/
#   HomoSapiens/
#    human_g1k_v37.fasta [2]
#
# REFERENCES:
# [1] http://aroma-project.org/data/annotationData/chipTypes/GenericHuman/
# [2] ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/technical/reference/
#                                                   human_g1k_v37.fasta.gz
############################################################################
library("aroma.seq")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
fullTest <- fullTest && isDirectory("annotationData,aroma.seq,private")
fullTest <- fullTest && isDirectory("fastqData,aroma.seq,private")
fullTest <- fullTest && isPackageInstalled("BatchJobs")
if (fullTest) {

setOption(aromaSettings, "devel/parallel", "BiocParallel")

dataSet <- "AlbertsonD_2012-SCC,AB042";
organism <- "HomoSapiens";

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism, prefix="human_g1k_v37")
print(fa)

# Binning reference (by 50kb)
ugp <- AromaUgpFile$byChipType("GenericHuman", tags="50kb")
print(ugp)

# Binned nucleotide-content reference
# (created from BSgenome.Hsapiens.UCSC.hg19 Bioconductor reference)
unc <- getAromaUncFile(ugp)
print(unc)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("fastqData,aroma.seq,private", dataSet, organism)
fqs <- IlluminaFastqDataSet$byPath(path)
fqs <- extract(fqs, 1:2)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build BWA index set (iff missing)
is <- buildBwaIndexSet(fa, method="is", verbose=verbose)
print(is)

# In addition to SAM read group data inferred from the Illumina FASTQ
# files, manual set the library information for the whole data set.
setSamReadGroup(fqs, SamReadGroup(LB="MPS-034"))

# BWA with BWA 'aln' options '-n 2' and '-q 40'.
alg <- BwaAlignment(fqs, indexSet=is, n=2, q=40)
print(alg)

bams <- process(alg, verbose=verbose)
print(bams)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Count reads per bin
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bc <- TotalCnBinnedCounting(bams, targetUgp=ugp)
print(bc)

ds <- process(bc, verbose=verbose)
verbose && print(verbose, ds)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Normalize GC content (and rescale to median=2 assuming diploid)
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bgn <- BinnedGcNormalization(ds)
verbose && print(verbose, bgn)

dsG <- process(bgn, verbose=verbose)
verbose && print(verbose, dsG)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Segmentation of tumors and normals independently (without a reference)
# and generation of a Chromosome Explorer report
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
seg <- CbsModel(dsG, ref="constant(2)", maxNAFraction=2/3)
verbose && print(verbose, seg)

ce <- ChromosomeExplorer(seg)
verbose && print(verbose, ce)
process(ce, chromsomes=c(1,19), verbose=verbose)

} # if (fullTest)


############################################################################
# HISTORY:
# 2013-11-02
# o Added to package's system tests.
#   Adopted from testScripts/complete/VUMC/91.FASTQ-to-CBS.R.
############################################################################
