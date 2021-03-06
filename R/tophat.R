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
#   \item{bowtieRefIndexPrefix}{A @character string specifying the Bowtie2 reference index prefix.}
#   \item{reads1}{(required) A @vector of FASTQ pathnames of reads.}
#   \item{reads2}{(optional; paired-end only) A @vector of FASTQ pathnames of mate reads.}
#   \item{gtf}{(optional) A GTF pathname.}
#   \item{mateInnerDist, mateStdDev}{(optional; paired-end only) The expected mean and standard
#         deviation of the inner distance between mate pairs.}
#   \item{optionsVec}{Vector of named options to pass to the executable.}
#   \item{...}{(Not used)}.
#   \item{outPath}{Directory where result files are written.}
#   \item{command}{The name of the executable.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \section{Support for compressed input files}{
#   TopHat (>= 1.3.0), which was released June 2011, handles FASTQ files that have been compressed
#   by gzip (or bzip2) [1]. If not supported, this method will give an informative error message
#   about it.
# }
#
# @author "HB,TT"
#
# \references{
#  [1] TopHat - A spliced read mapper for RNA-Seq, 2015.
#      \url{http://ccb.jhu.edu/software/tophat/}
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("tophat", "default", function(bowtieRefIndexPrefix, reads1=NULL, reads2=NULL, gtf=NULL, transcriptomeIndexPrefix=NULL, mateInnerDist=NULL, mateStdDev=NULL, optionsVec=NULL, ..., outPath="tophat/", command="tophat", verbose=FALSE) {
  # Make sure to evaluate registered onExit() statements
  on.exit(eval(onExit()))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'outPath'
  outPath <- Arguments$getWritablePath(outPath, mustNotExist=TRUE)

  # Argument 'bowtieRefIndexPrefix'
  # (and the existence of the corresponding directory)
  bowtieRefIndexPrefix <- Arguments$getCharacter(bowtieRefIndexPrefix, length=c(1L,1L))
  bowtieRefIndexPath <- dirname(bowtieRefIndexPrefix)
  bowtieRefIndexName <- basename(bowtieRefIndexPrefix)
  bowtieRefIndexPath <- Arguments$getReadablePath(bowtieRefIndexPath, absolutePath=TRUE)

  # Argument 'reads1'
  if (length(reads1) > 0L) {
    reads1 <- Arguments$getReadablePathnames(reads1, absolutePath=TRUE)
    assertNoDuplicated(reads1)
  }

  # Argument 'reads2'
  isPaired <- (length(reads2) > 0L)
  if (isPaired) {
    .stop_if_not(length(reads2) == length(reads1))
    reads2 <- Arguments$getReadablePathnames(reads2, absolutePath=TRUE)
    assertNoDuplicated(reads2)
  }

  # Argument 'gtf'
  if (!is.null(gtf)) {
    gtf <- Arguments$getReadablePathname(gtf, absolutePath=TRUE)
    if (isGzipped(gtf)) {
      throw("TopHat does not support gzipped GTF files: ", gtf)
    }
  }

  # Argument 'transcriptomeIndexPrefix'
  # (and the existence of the corresponding directory)
  if (!is.null(transcriptomeIndexPrefix)) {
    transcriptomeIndexPrefix <- Arguments$getCharacter(transcriptomeIndexPrefix, length=c(1L,1L))
    transcriptomeIndexPath <- dirname(transcriptomeIndexPrefix)
    transcriptomeIndexName <- basename(transcriptomeIndexPrefix)
    transcriptomeIndexPath <- Arguments$getReadablePath(transcriptomeIndexPath, absolutePath=TRUE)
  }


  # Argument 'mateInnerDist' & 'mateStdDev':
  if (!is.null(mateInnerDist)) {
    if (!isPaired) {
      throw("Argument 'mateInnerDist' can only be used with paired-end reads.")
    }
    mateInnerDist <- Arguments$getInteger(mateInnerDist, range=c(0,Inf))
  }
  if (!is.null(mateStdDev)) {
    if (!isPaired) {
      throw("Argument 'mateStdDev' can only be used with paired-end reads.")
    }
    mateStdDev <- Arguments$getInteger(mateStdDev, range=c(0,Inf))
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
     pushState(verbose)
     on.exit(popState(verbose), add=TRUE)
  }

  # Nothing to do if no reads and no transcript model
  if (is.null(reads1) && is.null(gtf)) {
    throw("Arguments 'reads1' and 'gtf' cannot both be NULL")
  }

  verbose && enter(verbose, "Running tophat()")
  if (length(reads1) > 0L) {
    verbose && cat(verbose, "R1 FASTQ files:")
    verbose && print(verbose, sapply(reads1, FUN=getRelativePath))
    if (isPaired) {
      verbose && cat(verbose, "R2 FASTQ files:")
      verbose && print(verbose, sapply(reads2, FUN=getRelativePath))
    }
  }
  verbose && cat(verbose, "Bowtie2 reference index prefix: ", bowtieRefIndexPrefix)
  verbose && cat(verbose, "Output directory: ", outPath)
  verbose && cat(verbose, "TopHat executable: ", command)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check whether gzipped FASTQ files are supported or not
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  gzipped <- any(regexpr("[.]gz$", c(reads1, reads2), ignore.case=TRUE) != -1L)
  if (gzipped) {
    verbose && cat(verbose, "Detected gzip'ed FASTQ files.")
    bin <- findTopHat(command=command)
    if (is.null(bin)) throw("TopHat executable not available.")
    verbose && str(verbose, bin)
    ver <- attr(bin, "version")
    verbose && cat(verbose, "TopHat version: ", ver)
    if (ver < "1.3.0") {
      throw("Detected gzip'ed FASTQ files, which is only supported by TopHat (>= 1.3.0): ", ver)
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate output atomically by writing to a temporary directory
  # that is renamed upon completion.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  outPathOrg <- outPath
  outPath <- sprintf("%s.TMP", outPath)
  outPath <- Arguments$getWritablePath(outPath, mustNotExist=TRUE)
  verbose && cat(verbose, "Temporary output directory: ", outPath)
  # At the end, assume failure, unless successful.
  outPathFinal <- sprintf("%s.ERROR", outPath)
  onExit({
    removeDirectory(outPathOrg, recursive=FALSE, mustExist=FALSE)
    file.rename(outPath, outPathFinal)
  })


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # While running, let the output directory be the working directory.
  # Make sure the working directory is restored when exiting.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opwd <- setwd(outPath)
  verbose && cat(verbose, "Original working directory: ", opwd)
  onExit({
    verbose && cat(verbose, "Resetting working directory: ", opwd)
    setwd(opwd)
  })


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Workaround the fact that tophat2 binary does not support commas
  # nor spaces in input pathnames, e.g. reference index files and
  # FASTQ files.
  #
  # NOTE: The current workaround only adjusts for commas in the path
  # names, not in the filenames.  If the filenames have commas,
  # the below assertion tests will throw an error.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # (1b) Inside the temporary output directory, setup a temporary
  #      input directory without commas
  inPath <- Arguments$getWritablePath("src")
  onExit({ removeDirectory(inPath, recursive=FALSE, mustExist=FALSE) })

  # (2a) Link to the bowtie2 index directory
  #      (such that tophat sees no commas)
  link <- file.path(inPath, "refIndex")
  bowtieRefIndexPath <- createLink(link=link, target=bowtieRefIndexPath)
  onExit({ removeDirectory(bowtieRefIndexPath) })
  link <- NULL # Not needed anymore
  bowtieRefIndexPrefix <- file.path(bowtieRefIndexPath, basename(bowtieRefIndexPrefix))
  bowtieRefIndexPrefix <- Arguments$getTuxedoOption(bowtieRefIndexPrefix)

  # (2b) Link to the GTF file
  if (!is.null(gtf)) {
    link <- file.path(inPath, basename(gtf))
    gtf <- createLink(link=link, target=gtf)
    onExit({ file.remove(gtf) })
    link <- NULL # Not needed anymore
    gtf <- Arguments$getTuxedoOption(gtf)
  }

  # (2c) Link to the tophat2 transcriptome index directory
  #      (such that tophat sees no commas)
  if (!is.null(transcriptomeIndexPrefix)) {
    link <- file.path(inPath, "transIndex")
    transcriptomeIndexPath <- createLink(link=link, target=transcriptomeIndexPath)
    onExit({ removeDirectory(transcriptomeIndexPath) })
    link <- NULL # Not needed anymore
    transcriptomeIndexPrefix <- file.path(transcriptomeIndexPath, basename(transcriptomeIndexPrefix))
    transcriptomeIndexPrefix <- Arguments$getTuxedoOption(transcriptomeIndexPrefix)
  }

  # (3a) Link to the FASTQ 'R1'
  #      (such that tophat sees no commas)
  if (length(reads1) > 0L) {
    reads1 <- sapply(reads1, FUN=function(pathname) {
      link <- file.path(inPath, basename(pathname))
      assertNoCommas(link)
      createLink(link=link, target=pathname)
    })
    onExit({ file.remove(reads1) })
    reads1 <- Arguments$getTuxedoOption(reads1)
  }

  # (3b) Link to the (optional) FASTQ 'R2'
  #      (such that tophat sees no commas)
  if (length(reads2) > 0L) {
    reads2 <- sapply(reads2, FUN=function(pathname) {
      link <- file.path(inPath, basename(pathname))
      assertNoCommas(link)
      createLink(link=link, target=pathname)
    })
    onExit({ file.remove(reads2) })
    reads2 <- Arguments$getTuxedoOption(reads2)
  }


  # When reaching this point, 'tophat' will not be able to tell whether
  # we are using the original files or the temporary workaround ones.
  # This is also reflected in the variable names, such that the below
  # code would look the same regardless whether we use a workaround or
  # not.  This means that if tophat one day will support commas in
  # directories and filenames, we can just delete this whole section
  # and everything will work out of the box. /HB 2013-10-31


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup tophat arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # The output directory should always to be the current directory
  opts <- c("-o"=".")

  # Append optional arguments
  if (!is.null(gtf)) opts <- c(opts, "-G"=shQuote(gtf))
  if (!is.null(transcriptomeIndexPrefix)) opts <- c(opts, "--transcriptome-index"=shQuote(transcriptomeIndexPrefix))
  if (!is.null(mateInnerDist)) opts <- c(opts, "--mate-inner-dist"=mateInnerDist)
  if (!is.null(mateStdDev)) opts <- c(opts, "--mate-std-dev"=mateStdDev)

  # Append user options
  opts <- c(opts, optionsVec)

  # At the very end:
  # (a) Append the bowtie2 reference index prefix
  opts <- c(opts, shQuote(bowtieRefIndexPrefix))

  # (b) Append the R1 FASTQ files
  if (length(reads1) > 0L) opts <- c(opts, shQuote(paste(reads1, collapse=",")))

  # (c) Paired-end analysis?  Then append the R2 FASTQ files
  if (length(reads2) > 0L) opts <- c(opts, shQuote(paste(reads2, collapse=",")))

  # Assert no duplicated options
  names <- names(opts)
  names <- names[nchar(names) > 0L]
  dups <- names[duplicated(names)]
  if (length(dups) > 0L) {
    throw("Duplicated options detected: ", paste(sQuote(dups), collapse=", "))
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call TopHat executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Calling systemTopHat()")
  args <- list(command=command, args=opts)
  verbose && cat(verbose, "Arguments:")
  verbose && print(verbose, args)

  args$verbose <- less(verbose, 10)
  res <- do.call(systemTopHat, args=args)
  status <- attr(res, "status"); if (is.null(status)) status <- 0L
  verbose && cat(verbose, "Results:")
  verbose && str(verbose, res)
  verbose && cat(verbose, "Status:")
  verbose && str(verbose, status)
  verbose && exit(verbose)

  # Successful?
  if (status == 0L) {
    # If we get this far, assume it was all successful.

    # Allow the temporary output path to be renamed to the
    # intended output path instead of the "error" one.
    outPathFinal <- outPathOrg
  }

  verbose && exit(verbose)

  res
}) # tophat()


setMethodS3("tophat1", "default", function(..., command="tophat") {
  tophat(..., command=command)
})

setMethodS3("tophat2", "default", function(..., command="tophat2") {
  tophat(..., command=command)
})
