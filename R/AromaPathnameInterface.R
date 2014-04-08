setConstructorS3("AromaPathnameInterface", function(...) {
  extend(Interface(...), "AromaPathnameInterface");
})

setMethodS3("directoryStructure", "AromaPathnameInterface", function(this, default="<rootpath>/<dataset>/<organism>/<sample>/", ...) {
  if (is.null(default)) default <- .findDefaultDirectoryStructure(this);
  NextMethod("directoryStructure", default=default);
})

setMethodS3("getOrganism", "AromaPathnameInterface", function(this, ...) {
  directoryItem(this, name="organism");
})



setConstructorS3("AromaSeqDataFile", function(...) {
  extend(AromaPathnameInterface(...), "AromaSeqDataFile");
})

setMethodS3("getDefaultFullName", "AromaSeqDataFile", function(this, ...) {
  value <- directoryItem(this, name="sample", mustExist=FALSE);
  if (is.null(value)) {
    value <- NextMethod("getDefaultFullName");
  } else {
    pattern <- getExtensionPattern(this);
    value <- gsub(pattern, "", value);
  }
  value;
})


setConstructorS3("AromaSeqDataFileSet", function(...) {
  extend(AromaPathnameInterface(...), "AromaSeqDataFileSet");
})

setMethodS3("getDefaultFullName", "AromaSeqDataFileSet", function(this, ...) {
  directoryItem(this, name="dataset");
})


############################################################################
# HISTORY:
# 2014-04-07
# o Added AromaPathnameInterface.
# o Added (temporary) AromaSeqData(File|FileSet).
############################################################################
