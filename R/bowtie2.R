###########################################################################/**
# @RdocDefault bowtie2
# @alias bowtie2_hb
#
# @title "Calls the Bowtie2 executable to align input reads"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{bowtieRefIndexPrefix}{bowtie2 reference index (partial pathname, minus the .1.bt2 suffix)}
#   \item{reads1}{Vector of fastq files to align, paired with reads2}
#   \item{reads2}{Vector of fastq files to align, paired with reads1}
#   \item{readsU}{Vector of fastq files to align (at least one of reads1 or readsU must be non-null}
#   \item{samFile}{Output file name}
#   \item{optionsVec}{Vector of named options (do not include names x, 1, 2, U, or S)}
#   \item{...}{...}
#   \item{commandName}{Name of executable}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
#
# @author "TT, HB"
#
# @keyword internal
#*/###########################################################################

## TODO:  This function has not been tested; the logic may not be complete; etc.

setMethodS3("bowtie2", "default", function(bowtieRefIndexPrefix=NULL, ##  ## Index filename prefix (i.e. minus trailing .X.bt2)
                                           reads1=NULL,  ## vector of pathnames, #1 mates
                                           reads2=NULL,  ## vector of pathnames, #2 mates
                                           readsU=NULL,  ## vector of pathnames, unpaired reads
                                           samFile=NULL,  ## SAM file for output
                                           optionsVec,
                                           ...,
                                           commandName='bowtie2',
                                           verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'bowtieRefIndexPrefix'
  # - check for bowtie2 reference index  ## TODO: ADD SUPPORT FOR BOWTIE1 INDICES
  if (!is.null(bowtieRefIndexPrefix)) {
    bowtieRefIndex1 <- paste(bowtieRefIndexPrefix, ".1.bt2", sep="")  ## (<<< assumes bowtie2)
    bowtieRefIndex1 <- Arguments$getReadablePathname(bowtieRefIndex1);
  } else {
    throw("Argument bowtieRefIndexPrefix is empty; supply (prefix of) bowtie reference index")
  }

  # Argument 'reads1', 'reads2', 'readsU'
  if ((is.null(reads1) && is.null(reads2)) &&
      is.null(readsU)) {
    throw("Arguments ('reads1' and 'reads2') and 'readsU' cannot all be empty; specify paired and/or unpaired reads")
  }
  # Argument 'reads1' and 'reads2'
  if (!is.null(reads1)) {
    reads1 <- sapply(reads1, FUN=Arguments$getReadablePathname)
    if (!is.null(reads2)) {
      reads2 <- sapply(reads2, FUN=Arguments$getReadablePathname)
    } else {
      throw("Argument 'reads2' is empty; supply reads2 when using reads1 (or just supply readsU)")
    }
  }
  # Argument 'readsU'
  if (!is.null(readsU)) {
    readsU <- sapply(readsU, FUN=Arguments$getReadablePathname)
  }

  # Argument 'samFile'
  if (!is.null(samFile)) {
    samFile <- Arguments$getWritablePathname(samFile)
  } else {
    throw("Argument 'samFile' is empty; supply an output file name")
  }

  ## Combine the above into "bowtie2 arguments"
  ## - Cf. usage:  "bowtie2 [options]* -x <bt2-idx> {-1 <m1> -2 <m2> | -U <r>} [-S <sam>]"
  bowtie2Args <- NULL  ## bowtie2 does not use 'arguments', just 'options'
  bowtie2Options <- c(x=bowtieRefIndexPrefix)
  if (!is.null(reads1)) {
    bowtie2Options <- c(bowtie2Options, '1'=unname(reads1))
    if (!is.null(reads2)) {
      bowtie2Options <- c(bowtie2Options, '2'=unname(reads2))
    }
  }
  if (!is.null(readsU)) {
    bowtie2Options <- c(bowtie2Options, 'U'=unname(reads1))
  }
  bowtie2Options <- c(bowtie2Options, 'S'=unname(samFile))

  ## Add dashes as appropriate to names of "bowtie2 options"
  bowtie2Options <- c(optionsVec, bowtie2Options)
  nms <- names(bowtie2Options)
  names(bowtie2Options) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")

  res <- do.call(what=systemBowtie2, args=list(command=commandName, args=c(bowtie2Options, bowtie2Args)))

  return(res)
})

############################################################################
# HISTORY:
# 2013-03-08
# o TT:  Completely rewritten to follow tophat template
# 2012-07-11
# o Created bowtie2() stub.
############################################################################
