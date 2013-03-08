###########################################################################/**
# @RdocDefault tophat
#
# @title "Calls the TopHat executable to align input reads"
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
setMethodS3("tophat", "default", function(command='tophat',
                                          bowtieRefIndexPrefix=NULL,  ## partial pathname (e.g. append .1.bt2 to get to a real file)
                                          reads1=NULL,  ## vector of pathnames
                                          reads2=NULL,  ## vector of pathnames
                                          optionsVec, ## vector of named options
                                          ..., verbose=FALSE) {

  ## ( Support a call like this: "tophat <options> bowtieRefIndexPrefix reads1 reads2" )

  # Argument 'bowtieRefIndexPrefix'
  # - check for bowtie2 reference index  ## TODO: ADD SUPPORT FOR BOWTIE1 INDICES
  if (!is.null(bowtieRefIndexPrefix)) {
    bowtieRefIndex1 <- paste(bowtieRefIndexPrefix, ".1.bt2", sep="")  ## (<<< assumes bowtie2)
    bowtieRefIndex1 <- Arguments$getReadablePathname(bowtieRefIndex1);
  } else {
    throw("Argument bowtieRefIndexPrefix is empty; supply (prefix of) bowtie reference index")
  }

  # Argument 'reads1'
  if (!is.null(reads1))
    {
      reads1 <- sapply(reads1, FUN=Arguments$getReadablePathname)
    } else {
      throw("Argument reads1 is empty; supply at least one input read file")
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
  tophatOptions <- optionsVec
  nms <- names(tophatOptions)
  names(tophatOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")

  res <- do.call(what=systemTopHat, args=list(command=command, args=c(tophatOptions, tophatArgs)))

  return(res)
})


############################################################################
# HISTORY:
# 2013-03-07
# o TT: Changed interface (standardized argument names, set NULL defaults)
# 2013-02-08
# o TT:  Created
############################################################################


