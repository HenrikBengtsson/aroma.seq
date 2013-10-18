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
#   \item{initialTopHatOutDir}{(optional) TopHat output dir used by tophat executable}
#   \item{finalTopHatOutDir}{(optional) TopHat output dir when this tophat() method exits}
#   \item{optionsVec}{Vector of named options to pass to tophat}
#   \item{...}{(Not used)}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("tophat", "default", function(bowtieRefIndexPrefix=NULL,
                                          reads1=NULL,
                                          reads2=NULL,
                                          tophatInDir="tophatIn",
                                          initialTopHatOutDir="tophatAroma",
                                          finalTopHatOutDir=NULL,
                                          optionsVec=NULL,
                                          ## systems2ArgsList
                                          ..., command='tophat',
                                          verbose=FALSE) {

  ## ( Support a call like this: "tophat <options> bowtieRefIndexPrefix reads1 reads2" )

  # Argument 'initialTopHatOutDir'
  initialTopHatOutDir <- Arguments$getWritablePathname(initialTopHatOutDir)

  # # Argument 'finalTopHatOutDir'
  if (!is.null(finalTopHatOutDir)) {
    finalTopHatOutDir <- Arguments$getWritablePathname(finalTopHatOutDir)
  }

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
      # Use symbolic link to avoid TopHat issues w/ commas; this will still fail if the basename has commas
      symLink <- createLink(link=file.path(tophatInDir, basename(reads1)), reads1)
      reads1ForTopHat <- symLink
    } else {
      throw("Argument reads1 is empty; supply at least one input read file")
    }

  # Argument 'reads2'
  if (!is.null(reads2))
    {
      reads2 <- sapply(reads2, FUN=Arguments$getReadablePathname)
      symLink <- createLink(link=file.path(tophatInDir, basename(reads2)), reads2)
      reads2ForTopHat <- symLink
    }

  ## Combine the above into "tophat arguments"
  if (!is.null(reads2ForTopHat)) {
    tophatArgs <- c(bowtieRefIndexPrefix,
                    paste(reads1ForTopHat, collapse=","),
                    paste(reads2ForTopHat, collapse=","))
  } else {  ## Technically just the above is probably ok
    tophatArgs <- c(bowtieRefIndexPrefix,
                    paste(reads1ForTopHat, collapse=","))
  }

  ## Add dashes as appropriate to names of "tophat options"
  tophatOptions <- NULL
  optionsVec <- c(optionsVec, "o"=initialTopHatOutDir)
  if (!is.null(optionsVec)) {
    tophatOptions <- optionsVec
    nms <- names(tophatOptions)
    names(tophatOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")
  }

  res <- do.call(what=systemTopHat, args=list(command=command, args=c(tophatOptions, tophatArgs)))

  if (!is.null(finalTopHatOutDir)) {
    # Move TopHat output to 'final' (e.g. aroma.seq-specified) directory.
    # Do some gymnastics to get list of files and subdirs w/o the parent dir name
    dir0 <- getwd()
    setwd(initialTopHatOutDir)
    tophatOutFiles <- findFiles(path=".", firstOnly=FALSE, recursive=TRUE)
    setwd(dir0)
    sapply(tophatOutFiles, function(f) {renameFile(file.path(initialTopHatOutDir, f),
                                                   file.path(finalTopHatOutDir, f), recursive=TRUE)})
    removeDirectory(initialTopHatOutDir, recursive=TRUE)

  }

  # Clean up symbolic links to input files
  removeDirectory(tophatInDir, recursive=TRUE)

  return(res)
})


############################################################################
# HISTORY:
# 2013-03-07
# o TT: Changed interface (standardized argument names, set NULL defaults)
# 2013-02-08
# o TT:  Created
############################################################################


