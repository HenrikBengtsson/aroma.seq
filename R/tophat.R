###########################################################################/**
# @RdocGeneric tophat
# @alias tophat1
# @alias tophat2
# @alias tophat.default
# @alias tophat1.default
# @alias tophat2.default
#
# @title "Calls the TopHat executable to align input reads"
#
# \description{
#  @get "title".
# }
#
# \usage{
#  @usage tophat,default
#  @usage tophat1,default
#  @usage tophat2,default
# }
#
# \arguments{
#   \item{command}{Name of executable}
#   \item{bowtieRefIndexPrefix}{bowtie2 reference index (partial pathname, i.e. minus the .x.bt2 suffix)}
#   \item{reads1}{(required) Vector of fastq filenames to align; currently only a single filename is supported}
#   \item{reads2}{(optional) Vector of fastq filenames to align, paired with reads1; currently only a single filename is supported}
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
setMethodS3("tophat", "default", function(bowtieRefIndexPrefix,
                                          reads1,
                                          reads2=NULL,
                                          outDir=NULL,
                                          optionsVec=NULL,
                                          ...,
                                          .initialTopHatOutDir=NULL,
                                          command='tophat',
                                          verbose=FALSE) {

  # ( Support a command line like this: "tophat <options> bowtieRefIndexPrefix reads1 reads2" )
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hasCommas <- function(pathnames, ...) {
    (regexpr(",", pathnames, fixed=TRUE) != -1L);
  } # hasCommas()

  assertNoCommas <- function(pathnames, ...) {
    stopifnot(!any(hasCommas(pathnames)));
  } # assertNoCommas()
		       

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument '.initialTopHatOutDir'
  if (is.null(.initialTopHatOutDir)) {
    if (!is.null(outDir)) {
      .initialTopHatOutDir <- gsub(",", "_", outDir)
    } else {
      .initialTopHatOutDir <- paste("tophatOut", getChecksum(reads1), sep="")
    }
  }
  .initialTopHatOutDir <- Arguments$getWritablePathname(.initialTopHatOutDir)

  # Argument 'outDir'
  if (is.null(outDir)) {
    outDir <- .initialTopHatOutDir
  }
  outDir <- Arguments$getWritablePathname(outDir)

  # Dir for TopHat input symbolic links
  tophatInDir <- Arguments$getWritablePath(file.path(.initialTopHatOutDir, "tophatIn"))

  # Set up the reference index and input read filenames as symbolic links, to avoid TopHat issues w/ commas
  # (This will still fail if the basenames have commas)

  # Argument 'bowtieRefIndexPrefix'
  bowtieRefIndexPrefix <- Arguments$getCharacter(bowtieRefIndexPrefix, length=c(1L,1L))
  bowtieRefIndexDir <- dirname(bowtieRefIndexPrefix)
  bowtieRefIndexDir <- Arguments$getReadablePathname(bowtieRefIndexDir);
  bowtieRefIndexDirForTopHat <- createLink(link=file.path(tophatInDir, "bowtieRefIndexDir"), bowtieRefIndexDir)
  bowtieRefIndexPrefixForTopHat <- file.path(bowtieRefIndexDirForTopHat, basename(bowtieRefIndexPrefix))

  # Argument 'reads1'
  reads1 <- Arguments$getReadablePathname(reads1)
  symLink <- createLink(link=file.path(tophatInDir, basename(reads1)), reads1)
  reads1ForTopHat <- symLink

  # Argument 'reads2'
  reads2ForTopHat <- NULL
  if (!is.null(reads2))
  {
    reads2 <- Arguments$getReadablePathname(reads2)
    symLink <- createLink(link=file.path(tophatInDir, basename(reads2)), reads2)
    reads2ForTopHat <- symLink
  }



  # Check that input files to tophat executable do not have commas
  assertNoCommas(bowtieRefIndexPrefixForTopHat)
  assertNoCommas(reads1ForTopHat)
  assertNoCommas(reads2ForTopHat)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call the tophat executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Combine the above into "tophat arguments"
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

  # Call tophat
  res <- do.call(what=systemTopHat, args=list(command=command, args=c(tophatOptions, tophatArgs)))


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Cleanup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Clean up input symbolic links
  removeDirectory(tophatInDir, recursive=TRUE)
  
  # Move TopHat output to 'final' (e.g. aroma.seq-specified) directory
  if (!is.null(outDir) && !(identical(outDir, .initialTopHatOutDir))) {
    file.rename(.initialTopHatOutDir, outDir)
  }

  res
}) # tophat()


setMethodS3("tophat1", "default", function(..., command="tophat") {
  tophat(..., command=command);
})

setMethodS3("tophat2", "default", function(..., command="tophat2") {
  tophat(..., command=command);
})


############################################################################
# HISTORY:
# 2013-10-30 [HB]
# o Added documentation for tophat1() and tophat2().
# 2013-10-19 [HB]
# o Added tophat1() and tophat2() which are wrappers for tophat().
# 2013-03-07
# o TT: Changed interface (standardized argument names, set NULL defaults)
# 2013-02-08
# o TT:  Created
############################################################################
