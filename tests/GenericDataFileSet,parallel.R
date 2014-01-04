library("R.filesets")
library("aroma.seq")

# - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up a file set
# - - - - - - - - - - - - - - - - - - - - - - - -
path <- system.file(package="R.filesets")
ds <- GenericDataFileSet$byPath(path)

# - - - - - - - - - - - - - - - - - - - - - - - -
# Get the size of each file
# - - - - - - - - - - - - - - - - - - - - - - - -
# Alt 1.
res1 <- lapply(ds, FUN=getFileSize)
print(res1)

# Alt 2. (FIXME: This returns NULL)
res2 <- dsApply(ds, FUN=getFileSize, .parallel="none")
print(res2)

# Alt 3.
if (isPackageInstalled("BatchJobs")) {
  res3 <- dsApply(ds, FUN=getFileSize, .parallel="BiocParallel::BatchJobs")
  names(res3) <- getNames(ds)
  print(res3)
  stopifnot(identical(res3, res1))
}
