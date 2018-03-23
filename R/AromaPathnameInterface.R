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
    # Inferred from a filename?
    if (!isTRUE(attr(value, "hasTail"))) {
      pattern <- getExtensionPattern(this);
      value <- gsub(pattern, "", value);
    }
  }
  value;
})

setMethodS3("getDefaultSamReadGroup", "AromaSeqDataFile", function(this, ...) {
  SamReadGroup()
})


setMethodS3("setSamReadGroup", "AromaSeqDataFile", function(this, rg, ...) {
  # Argument 'rg':
  if (!is.null(rg)) rg <- Arguments$getInstanceOf(rg, "SamReadGroup")
  this$.rg <- rg
  invisible(this)
})

setMethodS3("getSamReadGroup", "AromaSeqDataFile", function(this, ...) {
  rg <- this$.rg
  if (is.null(rg)) rg <- getDefaultSamReadGroup(this, ...)
  rg
})

setConstructorS3("AromaSeqDataFileSet", function(...) {
  extend(AromaPathnameInterface(...), "AromaSeqDataFileSet");
})

setMethodS3("getDefaultFullName", "AromaSeqDataFileSet", function(this, ...) {
  directoryItem(this, name="dataset");
})
