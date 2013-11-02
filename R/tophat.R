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
#   \item{outPath}{Directory where result files are written}
#   \item{optionsVec}{Vector of named options to pass to tophat}
#   \item{...}{(Not used)}
#   \item{command}{Name of executable}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \section{Support for compressed input files}{
#   TopHat (>= 1.3.0) handles FASTQ files that have been compressed by gzip (or bzip2) [1].
#   If not supported, this method will give an informative error message about it.
# }
#
# @author "HB,TT"
#
# \references{
#  [1] TopHat, University of Maryland, 2013.
#      \url{http://http://tophat.cbcb.umd.edu/}
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("tophat", "default", function(bowtieRefIndexPrefix,
                                          reads1,
                                          reads2=NULL,
                                          outPath="tophat/",
                                          optionsVec=NULL,
                                          ...,
                                          command="tophat",
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
  # Argument 'outPath'
  # Here we assert that the output directory does not already exists,
  # that it can be created, but we really don't want to create it because
  # later the temporary directory will be renamed to it, so we delete
  # immediately.
  outPath <- Arguments$getWritablePath(outPath, mustNotExist=TRUE)
  unlink(outPath)

  # Argument 'bowtieRefIndexPrefix'
  # (and the existence of the corresponding directory)
  bowtieRefIndexPrefix <- Arguments$getCharacter(bowtieRefIndexPrefix,
                                                           length=c(1L,1L))
  bowtieRefIndexPath <- dirname(bowtieRefIndexPrefix)
  bowtieRefIndexPath <- Arguments$getReadablePath(bowtieRefIndexPath)

  # Argument 'reads1'
  reads1 <- Arguments$getReadablePathname(reads1)

  # Argument 'reads2'
  if (!is.null(reads2)) {
    reads2 <- Arguments$getReadablePathname(reads2)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check whether gzipped FASTQ files are supported or not
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  gzipped <- any(regexpr("[.]gz$", c(reads1, reads2), ignore.case=TRUE) != -1L);
  if (gzipped) {
    ver <- attr(findTopHat2(), "version");
    if (ver <= "1.3.0") {
      throw("Detected gzip'ed FASTQ files, which is only supported by TopHat (>= 1.3.0): ", ver);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Workaround the fact that tophat2 binary does not support commas
  # in input pathnames, e.g. reference index files and FASTQ files.
  #
  # NOTE: The current workaround only adjusts for commas in the path
  # names, not in the filenames.  If the filenames have commas,
  # the below assertion tests will throw an error.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Record the orginal output directory
  outPathOrg <- outPath;

  # (1a) Use a temporary output directory without commas
  outPath <- gsub(",", "_", outPath);
  # Record which subdirectories already exists, so we can know which
  # are create here and hence should be deleted at the end.
  depth <- 0L;
  while (!is.null(dir <- getParent(outPath, depth=depth))) {
    if (isDirectory(dir)) break;
    depth <- depth + 1L;
  }
  if (depth > 0L) {
    outPathToDelete <- getParent(outPath, depth=depth-1L);
  } else {
    outPathToDelete <- NULL;
  }
  outPath <- Arguments$getWritablePath(outPath);
  assertNoCommas(outPath);


  # (1b) Inside the temporary output directory, setup a temporary
  #      input directory without commas
  inPath <- file.path(outPath, "src")
  inPath <- Arguments$getWritablePath(inPath)
  assertNoCommas(inPath)


  # (2a) Link to the bowtie2 index directory
  #      (such that tophat sees no commas)
  link <- file.path(inPath, "refIndex")
  bowtieRefIndexPath <- createLink(link=link, bowtieRefIndexPath)
  # When done, make sure to remove the temporary file link
  on.exit({
    if (file.exists(bowtieRefIndexPath)) {
      file.remove(bowtieRefIndexPath)
    }
  }, add=TRUE)
  bowtieRefIndexPrefix <- file.path(bowtieRefIndexPath, basename(bowtieRefIndexPrefix))
  assertNoCommas(bowtieRefIndexPrefix)

  # (3a) Link to the FASTQ 'R1'
  #      (such that tophat sees no commas)
  link <- file.path(inPath, basename(reads1))
  reads1 <- createLink(link=link, reads1)
  # When done, make sure to remove the temporary file link
  on.exit({
    if (file.exists(reads1)) {
      file.remove(reads1)
    }
  }, add=TRUE)
  assertNoCommas(reads1)

  # (3b) Link to the (optional) FASTQ 'R2'
  #      (such that tophat sees no commas)
  if (!is.null(reads2)) {
    reads2 <- Arguments$getReadablePathname(reads2)
    link <- file.path(inPath, basename(reads2))
    reads2 <- createLink(link=link, reads2)
    # When done, make sure to remove the temporary file link
    on.exit({
      if (file.exists(reads2)) {
        file.remove(reads2)
      }
    }, add=TRUE)
    assertNoCommas(reads2)
  }

  # (4a) Finally, when done, make sure to remove the temporary input
  #      directory which should be empty at this point ('recursive=FALSE')
  on.exit({
    removeDirectory(inPath, recursive=FALSE, mustExist=FALSE)
  }, add=TRUE)

  # (4b) ...and lastly, rename the temporary output directory to the
  #      final one.
  on.exit({
    if (!identical(outPathOrg, outPath)) {
      file.rename(outPath, outPathOrg)
    }
    if (!is.null(outPathToDelete)) {
      removeDirectory(outPathToDelete, recursive=TRUE, mustExist=FALSE)
    }
  }, add=TRUE)

  # When reaching this point, 'tophat' will not be able to tell whether
  # we are using the original files or the temporary workaround ones.
  # This is also reflected in the variable names, such that the below
  # code would look the same regardless whether we use a workaround or
  # not.  This means that if tophat one day will support commas in
  # directories and filenames, we can just delete this whole section
  # and everything will work out of the box. /HB 2013-10-31


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call the tophat executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opts <- optionsVec
  if (length(opts) > 0L) {
    # Add dashes as appropriate to names of "tophat options"
    nms <- names(opts)
    nms <- paste(ifelse(nchar(nms) == 1, "-", "--"), nms, sep="")
    names(opts) <- nms
  }

  # Append the output directory
  opts <- c(opts, "-o"=outPath)

  # Append the bowtie2 reference index prefix
  opts <- c(opts, bowtieRefIndexPrefix)

  # Append the R1 FASTQ files
  opts <- c(opts, paste(reads1, collapse=","))

  # Paired-end analysis?  Then append the R2 FASTQ files
  if (!is.null(reads2)) {
    opts <- c(opts, paste(reads2, collapse=","))
  }

  # Call TopHat executable
  args <- list(command=command, args=opts)
  res <- do.call(what=systemTopHat, args=args)


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
