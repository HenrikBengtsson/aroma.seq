library("aroma.seq")
setOption("R.filesets/parallel", NULL)

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

# Alt 2. (default, i.e. "none")
res2 <- dsApply(ds, FUN=getFileSize)
print(res2)
stopifnot(identical(res2, res1))

# Alt 3. (via an internal loop)
res2 <- dsApply(ds, FUN=getFileSize, .parallel="none")
print(res2)
stopifnot(identical(res2, res1))

# Alt 4. (via BatchJobs)
if (isPackageInstalled("BatchJobs")) {
  res3 <- dsApply(ds, FUN=getFileSize, .parallel="BiocParallel::BatchJobs")
  print(res3)
  stopifnot(identical(res3, res1))
}

# Alt 5. (via BiocParallel)
if (isPackageInstalled("BiocParallel")) {
  res4 <- dsApply(ds, FUN=getFileSize, .parallel="BiocParallel")
  print(res4)
  stopifnot(identical(res4, res1))
}
