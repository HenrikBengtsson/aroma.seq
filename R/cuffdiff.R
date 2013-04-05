###########################################################################/**
# @RdocDefault cuffdiff
#
# @title "Calls cuffdiff on input bam file(s)"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{transcriptsGtf}{Gtf file of transcripts ('gene model', e.g. from TopHat)}
#   \item{bams}{Vector of pathnames (.sam or .bam)}
#   \item{optionsVec}{Vector of named options}
#   \item{system2ArgsList}{Named list of arguments to pass to internal system2 call.}
#   \item{...}{...}
#   \item{verbose}{See @see "R.utils::Verbose".}
#   \item{commandName}{Name of executable}


  # }

#
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


  ## Build cuffdiff argument string
  ## - Cf. usage:   cuffdiff [options] <transcripts.gtf> <sample1_hits.sam> <sample2_hits.sam> [... sampleN_hits.sam]

  ## Construct cuffdiff 'arguments' (as opposed to 'options')
  cuffdiffArgs <- c(transcriptsGtf, as.vector(bams))

  if (FALSE) {  ## TAT: Superseded by more complex version below; do we need the extra checks?
    ## Add dashes as appropriate for cuffdiff options
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
        throw("Detected non-valid command line switch(es): ", hpaste(cuffdiffOptions[missing][invalid]));
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
