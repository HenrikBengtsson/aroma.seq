setConstructorS3("GcBaseFile", function(...) {
  extend(GenericDataFile(...), c("GcBaseFile", uses("AromaSeqDataFile")))
})


setMethodS3("getSeqNames", "GcBaseFile", function(this, ..., force=FALSE) {
  key <- list(method="getSeqNames", class=class(this), getPathname(this))
  dirs <- c("aroma.seq", "annotationData")
  res <- loadCache(key=key, dirs=dirs)
  if (!force && !is.null(res)) {
    return(res)
  }
  
  ## OS X requires 'zcat < infile'. /HB 2015-04-07
  bin <- if (isGzipped(this)) "zcat <" else "cat"
  cmd <- sprintf('%s "%s" | grep variableStep | cut -d " " -f 2', bin, getPathname(this))
  chrs <- system(cmd, intern=TRUE)
  chrs <- gsub(".*=", "", chrs)
  
  saveCache(chrs, key=key, dirs=dirs)
  chrs
})

setMethodS3("isCompatibleWith", "GcBaseFile", function(this, other, ...) {
  res <- NextMethod("isCompatibleWith")
  if (!isTRUE(res)) return(res)
  ## Assert same order
  idxs <- match(getSeqNames(other), getSeqNames(this))
  res <- all(diff(idxs) > 0, na.rm=TRUE)
  if (!res) {
    attr(res, "reason") <- "The ordering of sequence names does not match."
  }
  res
})
