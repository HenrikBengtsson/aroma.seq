###########################################################################/**
# @RdocDefault cufflinks
#
# @title "Calls cufflinks on input bam file(s)"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
# }
#
#
# @author
#*/###########################################################################
setMethodS3("cufflinks", "default", function(bams=NULL,  ## vector of pathnames
                                             optionsVec, ## vector of named options
                                             ...,
                                             commandName='cufflinks',
                                             verbose=FALSE) {

  ## ( Support a call like this: "cufflinks <options> bams" )

  # Argument 'bams'
  bams <- sapply(bams, FUN=Arguments$getReadablePathname)
  cufflinksArgs <- as.vector(bams)  ##  (DEV: UNLIST BETTER?)

  ## Add dashes as appropriate to names of "cufflinks options"
  cufflinksOptions <- optionsVec
  nms <- names(optionsVec)
  names(cufflinksOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")

  res <- do.call(what=systemCufflinks, args=list(commandName=commandName, args=c(cufflinksOptions, cufflinksArgs)))

  return(res)
})



