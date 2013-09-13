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
#   \item{command}{Name of executable}
#   \item{samFile}{sam file containing aligned reads, sorted by *name* if paired end}
#   \item{gfFile}{gene feature file, in gff format}
#   \item{outFile}{name of file to receive htseq-count output}
#   \item{optionsVec}{Vector of named options to pass to htseq-count}
#   \item{...}{(Not used)}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("htseqCount", "default", function(samFile=NULL,
                                              gfFile=NULL,
                                              outFile=NULL,
                                              optionsVec=c(s='no', a='10'),  ## Anders et al default
                                              ...,
                                              command='htseq-count',
                                              verbose=FALSE) {

  ## ( Support a call like this: "htseq-count -s no -a 10 lib_sn.sam gfFile > countFile")

  # Argument 'samFile'
  if (!is.null(samFile)) {
    samFile <- Arguments$getReadablePathname(samFile);
  } else {
    throw("Argument samFile is empty; supply sam file, sorted by name")
  }

  # Argument 'gfFile'
  if (!is.null(gfFile)) {
    gfFile <- Arguments$getReadablePathname(gfFile);
  } else {
    throw("Argument gfFile is empty; supply gene ('transcript') model file")
  }

  # Argument 'outFile'
  if (!is.null(outFile)) {
    outFile <- Arguments$getWritablePathname(outFile);
  } else {
    throw("Argument outFile is empty; supply output file name\n")
  }


  htseqArgs <- c(samFile, gfFile)
  htseqArgs <- c(htseqArgs, " > ", outFile)

  ## Add dashes as appropriate to names of "htseq options"
  htseqOptions <- NULL
  if (!is.null(optionsVec)) {
    htseqOptions <- optionsVec
    nms <- names(htseqOptions)
    names(htseqOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")
  }

  res <- do.call(what=systemHTSeqCount, args=list(command=command, args=c(htseqOptions, htseqArgs)))

  return(res)
})


############################################################################
# HISTORY:
# 2013-05-31
# o TT:  Created
############################################################################
