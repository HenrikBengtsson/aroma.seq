path <- system.file("testScripts/R", package="aroma.seq")
pathname <- file.path(path, "downloadUtils.R")
source(pathname)

library("aroma.seq")
verbose <- Arguments$getVerbose(-8, timestamp=TRUE)

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Downloading annotation data")

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Annotation data
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "annotationData/organisms/Lambda_phage/"
filename <- "lambda_virus.fa"
pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE)

if (!isFile(pathname)) {
  downloadBowtie2ExampleData()
}

fa <- FastaReferenceFile(filename, path=path)
print(fa)

verbose && exit(verbose)

