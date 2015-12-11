setConstructorS3("MPileupFile", function(...) {
  extend(GenericDataFile(...), "MPileupFile")
})

setMethodS3("getDefaultFullName", "MPileupFile", function(this, ...) {
  fullname <- NextMethod("getDefaultFullName")

  ## Support gzipped files
  ## FIXME: Can't this be done by TabularTextFile? /HB 2015-12-07
  if (isGzipped(this)) {
    pattern <- getExtensionPattern(this)
    fullname <- gsub(pattern, "", fullname)
  }

  fullname
}, protected=TRUE)


setConstructorS3("MPileupFileSet", function(...) {
  extend(GenericDataFileSet(...), "MPileupFileSet")
})

setMethodS3("findByName", "MPileupFileSet", function(static, ..., pattern="[.]mpileup(|[.]gz)$", paths=c("mpileupData")) {
  NextMethod("findByName", paths=paths, pattern=pattern)
}, static=TRUE)
