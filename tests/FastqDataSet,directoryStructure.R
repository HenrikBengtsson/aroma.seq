library("aroma.seq")

# Setup (writable) local data directory structure
setupExampleData()

dataset <- "TopHat-example"
organism <- "LambdaPhage"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup original FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs <- FastqDataSet$byName(dataset, organism=organism)

# Validate directory structure
for (ii in seq_along(fqs)) {
  fq <- fqs[[ii]]
  print(directoryStructure(fq))
  items <- directoryItems(fq)
  print(items)
  stopifnot(directoryItem(fq, "dataset") == dataset)
  stopifnot(directoryItem(fq, "organism") == organism)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create FASTQ set with custom directory structure
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
dataset2 <- sprintf("%s,customPath", dataset);
path2 <- file.path("fastqData", dataset2, organism);
path2 <- Arguments$getWritablePath(path2);

for (ii in seq_along(fqs)) {
  fq <- fqs[[ii]]
  sample <- getFullName(fq)
  path2T <- file.path(path2, sample)
  path2T <- Arguments$getWritablePath(path2T)
  pathname2T <- file.path(path2T, getFilename(fq))
  createLink(pathname2T, target=getPathname(fq), skip=TRUE)
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup custom FASTQ set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
fqs2 <- FastqDataSet$byPath(path2, struct="<rootpath>/<dataset>/<organism>/<sample>/", recursive=TRUE)

# Validate directory structure
for (ii in seq_along(fqs2)) {
  fq <- fqs2[[ii]]
  print(directoryStructure(fq))
  items <- directoryItems(fq)
  print(items)
  stopifnot(gsub(",.*", "", directoryItem(fq, "dataset")) == dataset)
  stopifnot(directoryItem(fq, "organism") == organism)
}
