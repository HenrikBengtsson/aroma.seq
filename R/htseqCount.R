###########################################################################/**
# @RdocDefault htseqCount
#
# @title "Calls the htseq-count executable to count input reads on features"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathnameS}{source file containing aligned reads; .sam and .bam supported}
#   \item{pathnameD}{(Optional) destination file to receive htseq-count output; use pathnameS with .count extension when NULL}
#   \item{gff}{gene feature file, in gff format}
#   \item{optionsVec}{Vector of named options to htseq-count}
#   \item{...}{(Not used)}
#   \item{command}{Name of executable}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("htseqCount", "default", function(pathnameS,
                                              pathnameD=NULL,
                                              gff,
                                              optionsVec=c('-s'='no', '-a'='10'),  ## Anders et al default
                                              ...,
                                              command='htseq-count',
                                              verbose=FALSE) {
  
  ## ( Support a call like this: "htseq-count -s no -a 10 lib_sn.sam gff > countFile")

  R.utils::use("Rsamtools")

  # BACKWARD COMPATIBILITY: Add asSam(), iff missing.
  if (packageVersion("Rsamtools") < "1.15.14") {
    asSam <- function(file, destination, ...) {
      file <- Arguments$getReadablePathname(file)
      fileD <- sprintf("%s.sam", destination)
      fileD <- Arguments$getWritablePathname(fileD)
      samtoolsView(file, fileD)
      fileD
    } # asSam()
  }

  # Argument 'pathnameS'
  pathnameS <- Arguments$getReadablePathname(pathnameS);
  
  # Argument 'pathnameD'
  if (is.null(pathnameD)) {
    pathnameD <- sub("[.](bam|sam)$", ".count", pathnameS, ignore.case=TRUE)
    if (pathnameD==pathnameS) { # should not happen
      pathnameD <- paste(pathnameS, ".count", sep="")
    }
  }
  stopifnot(pathnameD!=pathnameS)
  pathnameD <- Arguments$getWritablePathname(pathnameD, mustNotExist=TRUE);

  # Argument 'gff'
  gff <- Arguments$getReadablePathname(gff);
  
  ########
  # This methods uses somewhat convoluted sam / bam manipulations:
  # (sam input > ) bam input > sorted bam > sorted sam > htseq-count
  ########
  
  # Create/setup bam input for sorting
  if (regexpr(".*[.]sam$", pathnameS, ignore.case=TRUE) != -1) {  # match .sam
    inBam <- sub("[.]sam$", ".InTmp.bam", pathnameS, ignore.case=TRUE)
    inBam <- Arguments$getWritablePathname(inBam, mustNotExist=TRUE)
    asSam(pathnameS, inBam)
    on.exit({
      file.remove(inBam)
    }, add=TRUE)
  } else if (regexpr(".*[.]bam$", pathnameS, ignore.case=TRUE) != -1) {  # match .bam
    inBam <- pathnameS
  } else {
    throw("Not a known (.bam/.sam) file extension: ", pathnameS);
  }
  
  # Sort input bam by name
  inBamS <- sub("[.]bam$", ",ByNameTmp.bam", inBam)
  inBamS <- Arguments$getWritablePathname(inBamS, mustNotExist=TRUE)
  sortBam(inBam, sub("[.]bam$", "", inBamS), byQname=TRUE)
  
  # Convert sorted bam to sam (needed by htseq-count)
  inSamS <- sub("[.]bam$", ".sam", inBamS)
  inSamS <- Arguments$getWritablePathname(inSamS)
  samtoolsView(inBamS, inSamS)
  on.exit({
    file.remove(inBamS, inSamS)
  }, add=TRUE)
  
  htseqArgs <- c(inSamS, gff)
  htseqArgs <- c(htseqArgs, " > ", pathnameD)
  htseqOptions <- NULL
  if (!is.null(optionsVec)) {
    htseqOptions <- optionsVec
    # Deprecate this (dashes now required to be added by the user)
    # nms <- names(htseqOptions); names(htseqOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")
  }
  
  res <- do.call(what=systemHTSeqCount, args=list(command=command, args=c(htseqOptions, htseqArgs)))
  
  return(res)
})


############################################################################
# HISTORY:
# 2013-12-16
# o TT:  local asSam() added (needed for Rsamtools < 1.15.14)
# 2013-11-29
# o TT:  Internal sort-by-name added
# 2013-05-31
# o TT:  Created
############################################################################
