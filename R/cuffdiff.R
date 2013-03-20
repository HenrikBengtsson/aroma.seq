## Currently this is still dev /debug version

###########################################################################/**
# @Rdocdefault cuffdiff
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
setMethodS3("cuffdiff", "default", function(transcriptsGtf, ## e.g. from TopHat?
                                            bams=NULL,  ## vector of pathnames (will assume bam or sam works here)
                                            optionsVec=NULL, ## vector of named options
                                            ...,
                                            system2ArgsList=NULL,
                                            commandName='cuffdiff',
                                            verbose=FALSE) {

  ## ( Support a call like this: "cuffdiff <options> bams" )

  # Argument 'transcriptsGtf'
  transcriptsGtf <- Arguments$getReadablePathname(transcriptsGtf)

  # Argument 'bams'
  bams <- sapply(bams, FUN=Arguments$getReadablePathname)
  if (length(bams) < 2) {
    throw("'bams' argument too short; need at least two sam/bam files")
  }

  dotArgs <- list(...)
  ## - This will include arguments to systemCuffdiff itself, e.g .fake


  ## Construct cuffdiff 'arguments' (as opposed to 'options')
  cuffdiffArgs <- c(transcriptsGtf, as.vector(bams))

  ## Add dashes as appropriate for cuffdiff options
  if (FALSE) {
    cuffdiffOptions <- optionsVec
    nms <- names(cuffdiffOptions)
    names(cuffdiffOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")
  }

  ## Assign cuffdiffOptions; add dashes to names if needed
  if (!is.null(optionsVec)) {
    cuffdiffOptions <- optionsVec
    nms <- names(cuffdiffOptions);
    missing <- grep("^[^-]", nms);
    if (length(missing) > 0L) {
      invalid <- (nchar(cuffdiffOptions[missing]) < 1);
      if (any(invalid)) {
        throw("Detected non-valid command line switched: ", hpaste(cuffdiffOptions[missing][invalid]));
      }
      ## nms[missing] <- sprintf("-%s", nms[missing]);
      nmsMissing <- nms[missing]
      nmsMissingWithDash <- paste(ifelse(nchar(nmsMissing) == 1, "-", "--"), nmsMissing, sep="")
      nms[missing] <- nmsMissingWithDash
      names(cuffdiffOptions) <- nms
    }
  }
  cuffdiffAllArgs <- c(cuffdiffOptions, cuffdiffArgs)
  cuffdiffStr <- paste(names(cuffdiffAllArgs), cuffdiffAllArgs, sep=" ")

  ArgsList <- list(commandName=commandName, cuffdiffStr=cuffdiffStr)
  if (!is.null(system2ArgsList)) {
    ArgsList <- c(ArgsList, list(system2ArgsList=system2ArgsList))
  }
  ArgsList <- c(ArgsList, dotArgs)

  res <- do.call(what=systemCuffdiff, args=ArgsList)
  return(res)
})



