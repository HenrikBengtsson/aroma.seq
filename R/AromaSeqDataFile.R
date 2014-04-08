setConstructorS3("AromaSeqDataFile", function(...) {
  extend(GenericDataFile(...), "AromaSeqDataFile");
})

setConstructorS3("AromaSeqDataFileSet", function(...) {
  extend(GenericDataFileSet(...), "AromaSeqDataFileSet");
})

############################################################################
# HISTORY:
# 2014-04-07
# o Added (temporary) AromaSeqData(File|FileSet).
############################################################################
