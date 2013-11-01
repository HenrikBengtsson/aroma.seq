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
#   \item{bowtieRefIndexPrefix}{bowtie2 reference index (partial pathname, i.e. minus the .x.bt2 suffix)}
#   \item{reads1}{(required) Vector of fastq filenames to align; currently only a single filename is supported}
#   \item{reads2}{(optional) Vector of fastq filenames to align, paired with reads1; currently only a single filename is supported}
#   \item{outDir}{Directory where result files are written}
#   \item{optionsVec}{Vector of named options to pass to tophat}
#   \item{...}{(Not used)}
#   \item{command}{Name of executable}
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
                                          outDir="tophat/",
                                          optionsVec=NULL,
                                          ...,
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
  # Argument 'outDir'
  # Here we assert that the output directory does not already exists,
  # that it can be created, but we really don't want to create it because
  # later the temporary directory will be renamed to it, so we delete
  # immediately.
  outDir <- Arguments$getWritablePath(outDir, mustNotExist=TRUE)
  unlink(outDir)

  # Argument 'bowtieRefIndexPrefix'
  # (and the existence of the corresponding directory)
  bowtieRefIndexPrefix <- Arguments$getCharacter(bowtieRefIndexPrefix,
                                                           length=c(1L,1L))
  bowtieRefIndexDir <- dirname(bowtieRefIndexPrefix)
  bowtieRefIndexDir <- Arguments$getReadablePath(bowtieRefIndexDir)

  # Argument 'reads1'
  reads1 <- Arguments$getReadablePathname(reads1)

  # Argument 'reads2'
  if (!is.null(reads2)) {
    reads2 <- Arguments$getReadablePathname(reads2)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Workaround the fact that tophat2 binary does not support commas
  # in input pathnames, e.g. reference index files and FASTQ files.
  #
  # NOTE: The current workaround only adjusts for commas in the path
  # names, not in the filenames.  If the filenames have commas,
  # the below assertion tests will throw an error.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (1a) Use a temporary output directory without commas
  tophatOutDir <- gsub(",", "_", outDir)
  tophatOutDir <- Arguments$getWritablePath(tophatOutDir)
  assertNoCommas(tophatOutDir)
  # When done, make sure to rename it to the final one.
  on.exit({
    if (!identical(outDir, tophatOutDir)) {
      file.rename(tophatOutDir, outDir)
    }
  }, add=TRUE)


  # (1b) Inside the temporary output directory, setup a temporary
  #      input directory without commas
  tophatInDir <- file.path(tophatOutDir, "tophatIn")
  tophatInDir <- Arguments$getWritablePath(tophatInDir)
  assertNoCommas(tophatInDir)


  # (2a) Link to the bowtie2 index directory
  #      (such that tophat sees no commas)
  link <- file.path(tophatInDir, "bowtieRefIndexDir")
  bowtieRefIndexDirForTopHat <- createLink(link=link, bowtieRefIndexDir)
  # When done, make sure to remove the temporary file link
  on.exit({
    file.remove(bowtieRefIndexDirForTopHat)
  }, add=TRUE)
  bowtieRefIndexPrefixForTopHat <- file.path(bowtieRefIndexDirForTopHat, basename(bowtieRefIndexPrefix))
  assertNoCommas(bowtieRefIndexPrefixForTopHat)

  # (3a) Link to the FASTQ 'R1'
  #      (such that tophat sees no commas)
  link <- file.path(tophatInDir, basename(reads1))
  reads1ForTopHat <- createLink(link=link, reads1)
  # When done, make sure to remove the temporary file link
  on.exit({
    file.remove(reads1ForTopHat)
  }, add=TRUE)
  assertNoCommas(reads1ForTopHat)

  # (3b) Link to the (optional) FASTQ 'R2'
  #      (such that tophat sees no commas)
  reads2ForTopHat <- NULL
  if (!is.null(reads2)) {
    reads2 <- Arguments$getReadablePathname(reads2)
    link <- file.path(tophatInDir, basename(reads2))
    reads2ForTopHat <- createLink(link=link, reads2)
    # When done, make sure to remove the temporary file link
    on.exit({
      file.remove(reads2ForTopHat)
    }, add=TRUE)
    assertNoCommas(reads2ForTopHat)
  }

  # (4) Finally, when done, make sure to remove the temporary input
  #     directory which should be empty at this point ('recursive=FALSE')
  on.exit({
    removeDirectory(tophatInDir, recursive=FALSE)
  }, add=TRUE)


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
  optionsVec <- c(optionsVec, "o"=tophatOutDir)
  if (!is.null(optionsVec)) {
    tophatOptions <- optionsVec
    nms <- names(tophatOptions)
    names(tophatOptions) <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")
  }

  # Call tophat
  res <- do.call(what=systemTopHat, args=list(command=command, args=c(tophatOptions, tophatArgs)))


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
# 2013-10-31 [HB]
# o ROBUSTNESS: Now tophat() does a better job in making sure all
#   temporary directories and file links are removed regardless how
#   the function exits.
# o CLEANUP: Separated the validation of the arguments and the
#   workaround for dealing with commas.
# o CLEANUP: Dropped unnecessary argument '.initialTopHatOutDir'.
# o Now arguments in help appear in the same order as in the code.
# o Made arguments 'bowtieRefIndexPrefix' and 'reads1' mandatory.
#   Argument 'outDir' now defaults to 'tophat/' (and not NULL).
# 2013-10-30 [HB]
# o Added documentation for tophat1() and tophat2().
# 2013-10-19 [HB]
# o Added tophat1() and tophat2() which are wrappers for tophat().
# 2013-03-07
# o TT: Changed interface (standardized argument names, set NULL defaults)
# 2013-02-08
# o TT:  Created
############################################################################
