###########################################################################/**
# @RdocDefault cuffdiff
#
# @title "Calls cuffdiff on input bam file(s)"
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
setMethodS3("cuffdiff", "default", function(commandName='cuffdiff',
                                          bams=NULL,  ## vector of pathnames
                                          cuffdiffOptions, ## vector of named options
                                          ..., verbose=FALSE) {

  ## ( Support a call like this: "cuffdiff <options> bams" )

  # Argument 'bams'
  bams <- sapply(bams, FUN=Arguments$getReadablePathname)
  cuffdiffArgs <- as.vector(bams)  ##  (DEV: UNLIST BETTER?)

  ## Add dashes as appropriate to names of "cuffdiff options"
  nms <- names(cuffdiffOptions)
  names(cuffdiffOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")

  res <- do.call(what=systemCuffdiff, args=list(commandName=commandName, args=c(cuffdiffOptions, cuffdiffArgs)))

  return(res)
})



