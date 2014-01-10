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
setMethodS3("tophat", "default", function(bowtieRefIndexPrefix, reads1, reads2=NULL, gtf=NULL, outPath="tophat/", optionsVec=NULL, ..., command="tophat", verbose=FALSE) {
  # Make sure to evaluate registered onExit() statements
  on.exit(eval(onExit()));

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
  outPath <- Arguments$getWritablePath(outPath, mustNotExist=TRUE);

  # Argument 'bowtieRefIndexPrefix'
  # (and the existence of the corresponding directory)
  bowtieRefIndexPrefix <- Arguments$getCharacter(bowtieRefIndexPrefix,
                                                          length=c(1L,1L));
  bowtieRefIndexPath <- dirname(bowtieRefIndexPrefix);
  bowtieRefIndexName <- basename(bowtieRefIndexPrefix);
  bowtieRefIndexPath <- Arguments$getReadablePath(bowtieRefIndexPath, absolute=TRUE);

  # Argument 'reads1'
  stopifnot(length(reads1) > 0L);
  reads1 <- sapply(reads1, FUN=Arguments$getReadablePathname, absolute=TRUE);

  # Argument 'reads2'
  isPaired <- (length(reads2) > 0L);
  if (isPaired) {
    stopifnot(length(reads2) == length(reads1));
    reads2 <- sapply(reads2, FUN=Arguments$getReadablePathname, absolute=TRUE);
  }

  # Argument 'gtf'
  if (!is.null(gtf)) {
    gtf <- Arguments$getReadablePathname(gtf, absolute=TRUE);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
     pushState(verbose);
     on.exit(popState(verbose), add=TRUE);
  }


  verbose && enter(verbose, "Running tophat()");
  verbose && cat(verbose, "R1 FASTQ files:");
  verbose && print(verbose, sapply(reads1, FUN=getRelativePath));
  if (isPaired) {
    verbose && cat(verbose, "R2 FASTQ files:");
    verbose && print(verbose, sapply(reads2, FUN=getRelativePath));
  }
  verbose && cat(verbose, "Bowtie2 reference index prefix: ", bowtieRefIndexPrefix);
  verbose && cat(verbose, "Output directory: ", outPath);
  verbose && cat(verbose, "TopHat executable: ", command);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check whether gzipped FASTQ files are supported or not
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  gzipped <- any(regexpr("[.]gz$", c(reads1, reads2), ignore.case=TRUE) != -1L);
  if (gzipped) {
    verbose && cat(verbose, "Detected gzip'ed FASTQ files.");
    bin <- findTopHat(command);
    if (is.null(bin)) throw("TopHat executable not available.");
    verbose && str(verbose, bin);
    ver <- attr(bin, "version");
    verbose && cat(verbose, "TopHat version: ", ver);
    if (ver <= "1.3.0") {
      throw("Detected gzip'ed FASTQ files, which is only supported by TopHat (>= 1.3.0): ", ver);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate output atomically by writing to a temporary directory
  # that is renamed upon completion.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  outPathOrg <- outPath;
  outPath <- sprintf("%s.TMP", outPath);
  outPath <- Arguments$getWritablePath(outPath, mustNotExist=TRUE);
  verbose && cat(verbose, "Temporary output directory: ", outPath);
  # At the end, assume failure, unless successful.
  outPathFinal <- sprintf("%s.ERROR", outPath);
  onExit({
    removeDirectory(outPathOrg, recursive=FALSE, mustExist=FALSE);
    file.rename(outPath, outPathFinal);
  });


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # While running, let the output directory be the working directory.
  # Make sure the working directory is restored when exiting.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opwd <- setwd(outPath);
  verbose && cat(verbose, "Original working directory: ", opwd);
  onExit({
    verbose && cat(verbose, "Resetting working directory: ", opwd);
    setwd(opwd);
  });


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Workaround the fact that tophat2 binary does not support commas
  # in input pathnames, e.g. reference index files and FASTQ files.
  #
  # NOTE: The current workaround only adjusts for commas in the path
  # names, not in the filenames.  If the filenames have commas,
  # the below assertion tests will throw an error.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (1b) Inside the temporary output directory, setup a temporary
  #      input directory without commas
  inPath <- Arguments$getWritablePath("src");
  onExit({ removeDirectory(inPath, recursive=FALSE, mustExist=FALSE) })

  # (2a) Link to the bowtie2 index directory
  #      (such that tophat sees no commas)
  link <- file.path(inPath, "refIndex");
  bowtieRefIndexPath <- createLink(link=link, target=bowtieRefIndexPath);
  onExit({ file.remove(bowtieRefIndexPath) })
  bowtieRefIndexPrefix <- file.path(bowtieRefIndexPath, basename(bowtieRefIndexPrefix));
  assertNoCommas(bowtieRefIndexPrefix);

  # (2b) Link to the GTF file
  if (!is.null(gtf)) {
    link <- file.path(inPath, basename(gtf));
    gtf <- createLink(link=link, target=gtf);
    onExit({ file.remove(gtf) })
  }


  # (3a) Link to the FASTQ 'R1'
  #      (such that tophat sees no commas)
  reads1 <- sapply(reads1, FUN=function(pathname) {
    link <- file.path(inPath, basename(pathname));
    assertNoCommas(link);
    createLink(link=link, target=pathname);
  })
  onExit({ file.remove(reads1) })

  # (3b) Link to the (optional) FASTQ 'R2'
  #      (such that tophat sees no commas)
  if (isPaired) {
    reads2 <- sapply(reads2, FUN=function(pathname) {
      link <- file.path(inPath, basename(pathname));
      assertNoCommas(link);
      createLink(link=link, target=pathname);
    })
    onExit({ file.remove(reads2) })
  }


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

  # Set the output directory to be the current directory
  opts <- c(opts, "-o"=".")

  # GTF file, iff specified
  opts <- c(opts, "-G"=gtf)

  # Append the bowtie2 reference index prefix
  opts <- c(opts, bowtieRefIndexPrefix);

  # Append the R1 FASTQ files
  opts <- c(opts, paste(reads1, collapse=","));

  # Paired-end analysis?  Then append the R2 FASTQ files
  if (!is.null(reads2)) {
    opts <- c(opts, paste(reads2, collapse=","));
  }

  # Call TopHat executable
  verbose && enter(verbose, "Calling systemTopHat()");
  args <- list(command=command, args=opts);
  verbose && cat(verbose, "Arguments:");
  verbose && print(verbose, args);
  res <- do.call(systemTopHat, args=args);
  status <- attr(res, "status"); if (is.null(status)) status <- 0L;
  verbose && cat(verbose, "Results:");
  verbose && str(verbose, res);
  verbose && cat(verbose, "Status:");
  verbose && str(verbose, status);
  verbose && exit(verbose);

  # Successful?
  if (status == 0L) {
    # If we get this far, assume it was all successful.
    # Allow the temporary output path to be renamed to the
    # intended output path instead of the "error" one.
    outPathFinal <- outPathOrg;
  }

  verbose && exit(verbose);

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
# 2013-11-05 [HB]
# o Now the final output directory is <sampleName>.ERROR/, if the
#   run was not successful.
# 2013-11-02 [HB]
# o Utilizing onExit(), which is makes everything much easier.
# o Now the working directory is set to the output directory while
#   running TopHat.
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
