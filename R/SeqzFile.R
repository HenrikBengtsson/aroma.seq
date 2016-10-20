setConstructorS3("SeqzFile", function(..., .verify=FALSE) {
  extend(TabularTextFile(..., .verify=.verify), "SeqzFile")
})

setMethodS3("getDefaultFullName", "SeqzFile", function(this, ...) {
  fullname <- NextMethod("getDefaultFullName")

  ## Support gzipped files
  ## FIXME: Can't this be done by TabularTextFile? /HB 2015-12-07
  if (isGzipped(this)) {
    pattern <- getExtensionPattern(this)
    fullname <- gsub(pattern, "", fullname)
  }

  fullname
}, protected=TRUE)


setMethodS3("getDefaultColumnClassPatterns", "SeqzFile", function(this, ...) {
  c("*"="NULL", "(chromosome|base.ref|AB.|zygosity.normal|tumor.strand)"="character", "(position|depth.normal|depth.tumor|good.reads)"="integer", "(Af|Bf)|GC.percent"="double", "ratio"="NULL")
}, protected=TRUE)


setMethodS3("readAnnotationData", "SeqzFile", function(this, colClasses=c(chromosome="character", position="integer", GC.percent="numeric", base.ref="character"), ...) {
  readDataFrame(this, colClasses=colClasses, ...)
})


setConstructorS3("SeqzFileSet", function(...) {
  extend(TabularTextFileSet(...), "SeqzFileSet")
})

setMethodS3("findByName", "SeqzFileSet", function(static, ..., pattern="[.]seqz(|[.]gz)$", paths=c("seqzData")) {
  NextMethod("findByName", paths=paths, pattern=pattern)
}, static=TRUE)

setMethodS3("byName", "SeqzFileSet", function(static, ..., paths=c("seqzData")) {
  NextMethod("byName", paths=paths)
}, static=TRUE)

setMethodS3("readAnnotationData", "SeqzFileSet", function(this, colClasses=c(chromosome="character", position="integer", GC.percent="numeric", base.ref="character"), ...) {
  readDataFrame(this, colClasses=colClasses, ...)
})
