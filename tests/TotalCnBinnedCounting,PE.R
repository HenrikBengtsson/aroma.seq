library("aroma.seq")
setOption("R.filesets/parallel", "none")

fullTest <- (Sys.getenv("_R_CHECK_FULL_") != "")
fullTest <- fullTest && isCapableOf(aroma.seq, "bwa")
if (fullTest) {
# Setup (writable) local data directory structure
setupExampleData()

dataSet <- "YeastTest"
organism <- "SaccharomycesCerevisiae"

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTA reference file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fa <- FastaReferenceFile$byOrganism(organism)
print(fa)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataSet, organism=organism, paired=TRUE)
print(fqs)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Build Bowtie2 index set (iff missing)
is <- buildBowtie2IndexSet(fa, verbose=-10)
print(is)

# In addition to SAM read group data inferred from the Illumina FASTQ
# files, manual set the library information for the whole data set.
setSamReadGroup(fqs, SamReadGroup(LB="MPS-034"))

alg <- Bowtie2Alignment(fqs, indexSet=is)
print(alg)

bams <- process(alg, verbose=-10)
print(bams)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Generate target UGP on the fly
# NOTE: This is a rather ad hoc solution in order to be able
#       to tie into the aroma.cn framework. /HB 2013-11-12
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
chrLengths <- getSeqLengths(fa);
by <- 25e3
chrBins <- lapply(chrLengths, FUN=function(n) seq(from=1L, to=n, by=by))
nbrOfLoci <- sum(sapply(chrBins, FUN=length))

byTag <- sprintf("%gkb", by/1e3)
filename <- sprintf("%s.ugp", paste(c(organism, byTag), collapse=","))
ugp <- AromaUgpFile$allocate(filename=filename, path=tempfile(), nbrOfRows=nbrOfLoci, platform="HT-Seq", chipType=organism, overwrite=TRUE);
offset <- 0L
for (chr in seq_along(chrBins)) {
  x <- chrBins[[chr]]
  idxs <- offset + seq_along(x)
  ugp[idxs,1] <- chr
  ugp[idxs,2] <- x
  offset <- offset + length(x)
}
print(ugp)


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Count reads per bin
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bc <- TotalCnBinnedCounting(bams, targetUgp=ugp)
print(bc)

counts <- process(bc, verbose=-10)
verbose && print(verbose, counts)


data <- extractMatrix(counts, column=1L)
str(data)


} # if (fullTest)

############################################################################
# HISTORY:
# 2013-11-11
# o Created from BinnedGcNormalization test script.
############################################################################
