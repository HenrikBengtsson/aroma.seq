setConstructorS3("MPileupFile", function(...) {
  extend(GenericDataFile(...), "MPileupFile")
})

setConstructorS3("MPileupFileSet", function(...) {
  extend(GenericDataFileSet(...), "MPileupFileSet")
})
