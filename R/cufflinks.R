###########################################################################/**
# @RdocDefault cufflinks
#
# @title "Calls cufflinks on input bam file(s)"
#
# \description{
#  @get "title".
# }
#
# @synopsis  - level2
#
# \arguments{
# }
#
# \examples{\dontrun{
# }}
#
# @author
#*/###########################################################################
setMethodS3("cufflinks", "default", function(commandName='cufflinks',
                                          bams=NULL,  ## vector of pathnames
                                          cufflinksOptions, ## vector of named options
                                          ..., verbose=FALSE) {

  ## ( Support a call like this: "cufflinks <options> bams" )

  # Argument 'bams'
  bams <- sapply(bams, Arguments$getReadablePathname)
  cufflinksArgs <- as.vector(bams)

  ## Add dashes as appropriate to names of "cufflinks options"
  nms <- names(cufflinksOptions)
  names(cufflinksOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")

  res <- do.call(what="systemCufflinks", args=list(commandName=commandName, args=c(cufflinksOptions, cufflinksArgs)))

  return(res)
})



