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
#   \item{reads1}{(required) Vector of fastq filenames to align; currently only a single filename is supported}
#   \item{reads2}{(optional) Vector of fastq filenames to align, paired with reads1; currently only a single filename is supported}
#   \item{.tophatInDir}{(optional) Directory for symbolic links to input fastq files}
#   \item{.initialTopHatOutDir}{(optional) TopHat output dir used by tophat executable}
#   \item{outDir}{(optional) TopHat output dir at method exit}
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
                                          outDir=NULL,
                                          optionsVec=NULL,
                                          ...,
                                          .tophatInDir="tophatIn",
                                          .initialTopHatOutDir=NULL,
                                          command='tophat',
                                          verbose=FALSE) {
  
  # ( Support a call like this: "tophat <options> bowtieRefIndexPrefix reads1 reads2" )
  
  # Argument '.tophatInDir'
  .tophatInDir <- Arguments$getWritablePathname(.tophatInDir)
  
  # Argument '.initialTopHatOutDir'
  if (is.null(.initialTopHatOutDir)) {
    .initialTopHatOutDir <- gsub(",", "_", outDir)
  }
  .initialTopHatOutDir <- Arguments$getWritablePathname(.initialTopHatOutDir)

            
  # # Argument 'outDir'
  if (!is.null(outDir)) {
    outDir <- Arguments$getWritablePathname(outDir)
  }
  
  # Set up the reference index and input read filenames as symbolic links, to avoid TopHat issues w/ commas
  # (This will still fail if the basenames have commas)
  
  # Argument 'bowtieRefIndexPrefix'
  if (!is.null(bowtieRefIndexPrefix)) {
    bowtieRefIndexDir <- dirname(bowtieRefIndexPrefix)
    bowtieRefIndexDir <- Arguments$getReadablePathname(bowtieRefIndexDir);
    bowtieRefIndexDirForTopHat <- createLink(link=file.path(.tophatInDir, "bowtieRefIndexDir"), bowtieRefIndexDir)
    bowtieRefIndexPrefixForTopHat <- file.path(bowtieRefIndexDirForTopHat, basename(bowtieRefIndexPrefix))
  } else {
    throw("Argument bowtieRefIndexPrefix is empty; supply (prefix of) bowtie reference index")
  }

  reads1ForTopHat <- NULL
  reads2ForTopHat <- NULL
  # Argument 'reads1'
  if (!is.null(reads1))
  {
    reads1 <- Arguments$getReadablePathname(reads1)
    symLink <- createLink(link=file.path(.tophatInDir, basename(reads1)), reads1)
    reads1ForTopHat <- symLink
  } else {
    throw("Argument reads1 is empty; supply at least one input read file")
  }
  # Argument 'reads2'
  if (!is.null(reads2))
  {
    reads2 <- Arguments$getReadablePathname(reads2)
    symLink <- createLink(link=file.path(.tophatInDir, basename(reads2)), reads2)
    reads2ForTopHat <- symLink
  }
  
  # Check that input files to tophat executable do not have commas
  stopifnot(length(grep(",", c(bowtieRefIndexPrefixForTopHat, reads1ForTopHat, reads2ForTopHat))) == 0)
  
  ## Combine the above into "tophat arguments"
  if (!is.null(reads2ForTopHat)) {
    tophatArgs <- c(bowtieRefIndexPrefixForTopHat,
                    paste(reads1ForTopHat, collapse=","),
                    paste(reads2ForTopHat, collapse=","))
  } else {  ## Technically just the above is probably ok
    tophatArgs <- c(bowtieRefIndexPrefixForTopHat,
                    paste(reads1ForTopHat, collapse=","))
  }
  
  ## Add dashes as appropriate to names of "tophat options"
  tophatOptions <- NULL
  optionsVec <- c(optionsVec, "o"=.initialTopHatOutDir)
  if (!is.null(optionsVec)) {
    tophatOptions <- optionsVec
    nms <- names(tophatOptions)
    names(tophatOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")
  }
  
  res <- do.call(what=systemTopHat, args=list(command=command, args=c(tophatOptions, tophatArgs)))
    
  if (!is.null(outDir)) {
    # Move TopHat output to 'final' (e.g. aroma.seq-specified) directory.
    file.rename(.initialTopHatOutDir, outDir)
  }
  
  # Clean up symbolic links to input files
  removeDirectory(.tophatInDir, recursive=TRUE)
  
  # Remove temp dir
  # removeDirectory(dirname(.initialTopHatOutDir)) # Too dangerous as is; could do recursive delete of empty dirs
  
  return(res)
})


############################################################################
# HISTORY:
# 2013-03-07
# o TT: Changed interface (standardized argument names, set NULL defaults)
# 2013-02-08
# o TT:  Created
############################################################################


