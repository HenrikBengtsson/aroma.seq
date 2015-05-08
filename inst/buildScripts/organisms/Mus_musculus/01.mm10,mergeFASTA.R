# REFERENCES:
library("aroma.seq")

organism <- "Mus_musculus"
tag <- "mm10"
chrs <- c(1:19, "X", "Y", "M")
nchrs <- length(chrs)
filename <- sprintf("%s_%s_chr1-%d.fa", organism, tag, nchrs)

# Already done?
path <- file.path("annotationData", "organisms", organism)
pathname <- file.path(path, filename)
if (isFile(pathname)) {
  fa <- FastaReferenceFile(pathname)
} else {
  # Build!

  # Setup all FASTA reference files
  pathS <- file.path("annotationData", "organisms", organism, tag)
  mkdirs(pathS)
  fas <- FastaReferenceSet$byPath(pathS)

  # Extract the "base" ones
  fas <- extract(fas, "chr([0-9]{1,2}|[XYM])")

  # Download?
  if (length(fas) < 22L) {
    urlRoot <- "ftp://hgdownload.cse.ucsc.edu/goldenPath/mm10/chromosomes"
    filenames <- sprintf("chr%s.fa.gz", chrs)
    pathnames <- file.path(pathS, filenames)
    ok <- file_test("-f", pathnames)
    pathnames <- pathnames[!ok]
    for (kk in seq_along(pathnames)) {
      pathname <- pathnames[kk]
      url <- file.path(urlRoot, basename(pathname), fsep="/")
      message(url)
      pathname <- downloadFile(url, filename=pathname)
    }
    fas <- FastaReferenceSet$byPath(pathS)
    fas <- extract(fas, "chr([0-9]{1,2}|[XYM])")
  }
  print(fas)

  # Sanity check
  stopifnot(length(fas) == nchrs)

  # Order in which chromosomes are merged?
  # (a) chr1, ..., chr9, chr10, ..., chr19, chrM, chrX, chrY
  fas <- sortBy(fas, "mixedsort")

  # (b) chr1, chr10, chr11, ..., chr19, chr2, ..., chr9, chrM, chrX, chrY
  fas <- sortBy(fas, "lexicographic")

  # (c) chr10, chr11, ..., chr19, chr1, chr2, ..., chr9, chrM, chrX, chrY
  #     This is how it's done in ftp://ussd-ftp.illumina.com/Mus_musculus/
  #     UCSC/mm9/Mus_musculus_UCSC_mm9.tar.gz
  #     Doing it this way creates an identical file
  fas <- sortBy(fas, "lexicographic")
  fas <- extract(fas, order(nchar(getNames(fas)), decreasing=TRUE))

  # Merge FASTA files
  fa <- writeFastaReferenceFile(fas, filename=filename, path=path, verbose=TRUE)
}

print(fa)
faZ <- getChecksumFile(fa)
print(faZ)



############################################################################
# HISTORY:
# 2014-06-19 [HB]
# o Created and refined from mm9 version.
############################################################################
