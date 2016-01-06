###########################################################################/**
# @RdocClass BwaIndexSet
# @aliasmethod isCompatibleWith
# @aliasmethod isComplete
# @aliasmethod getSeqLengths
# @aliasmethod getSeqNames
# @aliasmethod readAnnData
# @aliasmethod as.character
#
# @title "The BwaIndexSet class"
#
# \description{
#  @classhierarchy
#
#  An BwaIndexSet object represents a set of @see "BwaIndexFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "AbstractIndexSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setConstructorS3("BwaIndexSet", function(...) {
  extend(AbstractIndexSet(...), c("BwaIndexSet", uses("SequenceContigsInterface")))
})


setMethodS3("as.character", "BwaIndexSet", function(x, ...) {
  s <- NextMethod("as.character")
  s <- c(s, getSeqGenericSummary(x, ...))
  s
})


setMethodS3("isComplete", "BwaIndexSet", function(this, ...) {
  knownExts <- c("amb", "ann", "bwt", "pac", "sa");
  if (length(this) < length(knownExts)) return(FALSE);

  exts <- sapply(this, getExtension);
  missing <- setdiff(knownExts, exts);
  if (any(missing)) {
    return(FALSE);
  }

  TRUE;
})


setMethodS3("getSeqLengths", "BwaIndexSet", function(this, ...) {
  data <- readAnnData(this)
  lengths <- data$length
  names(lengths) <- data$name
  lengths
})


setMethodS3("readAnnData", "BwaIndexSet", function(this, ...) {
  ## Identify sequence names and lengths
  pathnames <- getPathnames(this)
  ann <- is[[grep("[.]ann$", pathnames)]]

  pathname <- getPathname(ann)
  data <- readLines(pathname)

  ## Parse
  dataT <- lapply(data, FUN=function(x) {
    x <- gsub("^([^ ]+) ([^ ]+) (.+)", "\\1\t\\2\t\\3", x)
    strsplit(x, split="\t", fixed=TRUE)[[1]]
  })
  ns <- sapply(dataT, FUN=length)
  stopifnot(all(ns == 3))
  hdr <- dataT[[1]]
  hdr <- list(
    totalSeqLength=as.double(hdr[1]),
    nbrOfSequences=as.integer(hdr[2]),
    unknown=hdr[3]
  )
  dataT <- dataT[-1]
  stopifnot(length(dataT) == 2*hdr$nbrOfSequences)
  seqInfo <- list()
  for (kk in seq_len(hdr$nbrOfSequences)) {
    dataKK <- dataT[1:2 + 2*(kk-1)]
    dataKK <- unlist(dataKK, use.names=FALSE)
    seqInfoKK <- data.frame(a=dataKK[1], name=dataKK[2], description=dataKK[3], start=dataKK[4], length=as.double(dataKK[5]), c=dataKK[6], stringsAsFactors=FALSE)
    seqInfo[[kk]] <- seqInfoKK
  }
  seqInfo <- Reduce(rbind, seqInfo)
  stopifnot(sum(seqInfo$length) == hdr$totalSeqLength)
  attr(seqInfo, "header") <- hdr
  seqInfo
})



############################################################################
# HISTORY:
# 2012-09-25
# o Created.
############################################################################
