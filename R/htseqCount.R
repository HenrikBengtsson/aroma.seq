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
#   \item{inFile}{input file containing aligned reads, sorted by *name* if paired end; .sam and .bam supported}
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
setMethodS3("htseqCount", "default", function(inFile,
                                              gfFile,
                                              outFile=NULL,
                                              optionsVec=c(s='no', a='10'),  ## Anders et al default
                                              ...,
                                              command='htseq-count',
                                              verbose=FALSE) {

  ## ( Support a call like this: "htseq-count -s no -a 10 lib_sn.sam gfFile > countFile")

  # Argument 'inFile'
  inFile <- Arguments$getReadablePathname(inFile);

  # Assign 'samFile' = actual input filename, with support for .bam input
  if (regexpr(".*[.]sam$", inFile, ignore.case=TRUE) != -1) {
    samFile <- inFile
  } else if (regexpr(".*[.]bam$", inFile, ignore.case=TRUE) != -1) {
    samFile <- sub("[.]bam$", ".sam", inFile, ignore.case=TRUE)
    samFile <- Argument$getWritablePathname(samFile, mustNotExist=TRUE)
    samtoolsView(pathname=inFile, pathnameD=samFile)
    on.exit({
      file.remove(samFile)
    }, add=TRUE)
  } else {
    throw("Unknown file extension: ", inFile);
  }

  # Argument 'gfFile'
  gfFile <- Arguments$getReadablePathname(gfFile);

  # Argument 'outFile'
  outFile <- Arguments$getWritablePathname(outFile);

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
