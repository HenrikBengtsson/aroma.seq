###########################################################################/**
# @RdocDefault bowtie2Build
#
# @title "Creates index on reference genome using bowtie2-build"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{refReads}{Vector of (FASTA) files with reference sequence}
#   \item{bowtieRefIndexPrefix}{bowtie2 reference index to be built (partial pathname, minus the .*.bt2 suffix)}
#   \item{optionsVec}{Vector of named options}
#   \item{...}{...}
#   \item{commandName}{Name of executable}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
#
# @author
#*/###########################################################################
setMethodS3("bowtie2Build", "default", function(refReads=NULL,  ## vector of pathnames
                                                bowtieRefIndexPrefix=NULL, ## Index filename prefix (i.e. minus trailing .X.bt2)
                                                optionsVec=NULL, ## vector of named options
                                                ...,
                                                commandName='bowtie2-build',
                                                verbose=FALSE) {
  # Argument 'refReads'
  if (!is.null(refReads))
    {
      refReads <- sapply(refReads, FUN=Arguments$getReadablePathname)
    } else {
      throw("Argument 'refReads' is empty; supply at least one input read file")
    }

  # Argument 'bowtieRefIndexPrefix'
  # - check for bowtie2 reference index  ## TODO: ADD SUPPORT FOR BOWTIE1 INDICES
  if (!is.null(bowtieRefIndexPrefix)) {
    bowtieRefIndexPrefix <- Arguments$getWritablePathname(bowtieRefIndexPrefix);
  } else {
    throw("Argument bowtieRefIndexPrefix is empty; supply (prefix of) bowtie reference index")
  }

  ## Combine the above into "bowtie2 arguments"
  bowtie2Args <- c(paste(unname(refReads), collapse=","), bowtieRefIndexPrefix)   ## reference reads first, index base second

  ## Add dashes as appropriate to names of "bowtie2 options"
  bowtie2Options <- NULL
  if (!is.null(optionsVec)) {
    bowtie2Options <- optionsVec
    nms <- names(bowtie2Options)
    names(bowtie2Options) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")
  }

  res <- do.call(what=systemBowtie2Build, args=list(command=commandName, args=c(bowtie2Options, bowtie2Args)))

  ## DEV
  ## res <- do.call(what=systemBowtie2Build, args=list(command=commandName, args=c(bowtie2Options, bowtie2Args), verbose=TRUE, system2ArgsList=list(stderr=TRUE)))

  return(res)
})

############################################################################
# HISTORY:
# 2013-03-08
# o TT:  Completely rewritten to follow tophat.R template
# 2012-09-14
# o First real draft of "level 2" code (TAT)
# 2012-07-18
# o Created bowtie2-build() stub. (HB)
############################################################################
