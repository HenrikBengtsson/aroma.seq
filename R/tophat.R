###########################################################################/**
# @RdocDefault tophat
#
# @title "Calls tophat on input reads"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{...}
# }
#
# \examples{\dontrun{
# }}
#
# @author
#*/###########################################################################
setMethodS3("tophat", "default", function(commandName='tophat',
                                          bowtieRefIndexPrefix,  ## partial pathname (e.g. append .1.bt2 to get to a real file)
                                          reads1=NULL,  ## vector of pathnames
                                          reads2=NULL,  ## vector of pathnames
                                          tophatOptions, ## vector of named options
                                          ..., verbose=FALSE) {

  ## ( Support a call like this: "tophat <options> bowtieRefIndex reads1 reads2" )

  # Argument 'bowtieRefIndexPrefix'
  # - check for bowtie2 reference index  ## TODO: ADD SUPPORT FOR BOWTIE1 INDICES
  bowtieRefIndex1 <- paste(bowtieRefIndexPrefix, ".1.bt2", sep="")  ## (<<< assumes bowtie2)
  bowtieRefIndex1 <- Arguments$getReadablePathname(bowtieRefIndex1);

  # Argument 'reads1'
  if (!is.null(reads1))
    {
      reads1 <- sapply(reads1, FUN=Arguments$getReadablePathname)
    } else {
      throw("Argument reads1 is empty; need at least one input read")
    }

  # Argument 'reads2'
  if (!is.null(reads2))
    {
      reads2 <- sapply(reads2, FUN=Arguments$getReadablePathname)
    }

  ## Combine the above into "tophat arguments"
  if (!is.null(reads2)) {
    tophatArgs <- c(bowtieRefIndexPrefix, unname(reads1), unname(reads2))  ## (use unname to get rid of names)
  } else {  ## Technically just the above is probably ok
    tophatArgs <- c(bowtieRefIndexPrefix, unname(reads1))
  }

  ## Add dashes as appropriate to names of "tophat options"
  nms <- names(tophatOptions)
  names(tophatOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")

  res <- do.call(what=systemTophat, args=list(commandName=commandName, args=c(tophatOptions, tophatArgs)))

  return(res)
})



