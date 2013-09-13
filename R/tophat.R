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
#   \item{command}{Name of executable}
#   \item{bowtieRefIndexPrefix}{bowtie2 reference index (partial pathname, i.e. minus the .x.bt2 suffix)}
#   \item{reads1}{(required) Vector of fastq filenames to align}
#   \item{reads2}{(optional) Vector of fastq filenames to align, paired with reads1}
#   \item{optionsVec}{Vector of named options to pass to tophat}
#   \item{...}{(Not used)}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("tophat", "default", function(bowtieRefIndexPrefix=NULL,  ##
                                          reads1=NULL,
                                          reads2=NULL,
                                          optionsVec=NULL,
                                          ## systems2ArgsList
                                          ..., command='tophat',
                                          verbose=FALSE) {

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
    tophatArgs <- c(bowtieRefIndexPrefix,
                    paste(reads1, collapse=","),
                    paste(reads2, collapse=","))
  } else {  ## Technically just the above is probably ok
    tophatArgs <- c(bowtieRefIndexPrefix,
                    paste(reads1, collapse=","))
  }

  ## Add dashes as appropriate to names of "tophat options"
  tophatOptions <- NULL
  if (!is.null(optionsVec)) {
    tophatOptions <- optionsVec
    nms <- names(tophatOptions)
    names(tophatOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")
  }

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


