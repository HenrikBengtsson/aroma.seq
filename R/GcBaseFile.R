setConstructorS3("GcBaseFile", function(...) {
  extend(GenericDataFile(...), c("GcBaseFile", uses("AromaSeqDataFile"), uses("SequenceContigsInterface")))
})

setMethodS3("as.character", "GcBaseFile", function(x, ...) {
  s <- NextMethod("as.character")
  s <- c(s, getSeqGenericSummary(x, ...))
  s
})

setMethodS3("getSeqLengths", "GcBaseFile", function(this, ..., force=FALSE) {
  key <- list(method="getSeqLengths", class=class(this), getPathname(this))
  dirs <- c("aroma.seq", "annotationData")
  res <- loadCache(key=key, dirs=dirs)
  if (!force && !is.null(res)) {
    return(res)
  }

  ## OS X requires 'zcat < infile'. /HB 2015-04-07
  bin <- if (isGzipped(this)) "zcat <" else "cat"
  cmd <- sprintf('%s "%s" | grep variableStep | cut -d " " -f 2', bin, getPathname(this))
  names <- system(cmd, intern=TRUE)
  names <- gsub(".*=", "", names)

  lens <- rep(NA_integer_, times=length(names))
  names(lens) <- names

  saveCache(lens, key=key, dirs=dirs)

  lens
})
