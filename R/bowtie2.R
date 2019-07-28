###########################################################################/**
# @RdocFunction bowtie2
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
#   \item{reads1}{(required) A @vector of FASTQ pathnames of reads.}
#   \item{reads2}{(optional; paired-end only) A @vector of FASTQ pathnames of mate reads.}
#   \item{indexPrefix}{Bowtie2 reference index prefix.}
#   \item{pathnameSAM}{Output SAM file.}
#   \item{...}{...}
#   \item{gzAllowed}{A @logical specifying whether gzipped FASTQ files are supported or not.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns ...
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
bowtie2 <- function(reads1, reads2=NULL, indexPrefix, pathnameSAM, ..., gzAllowed=NA, verbose=FALSE) {
  # Make sure to evaluate registered onExit() statements
  on.exit(eval(onExit()))

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'reads1'
  reads1 <- Arguments$getReadablePathnames(reads1, absolutePath=TRUE)
  assertNoDuplicated(reads1)

  # Argument 'reads2'
  isPaired <- (length(reads2) > 0L)
  if (isPaired) {
    .stop_if_not(length(reads2) == length(reads1))
    reads2 <- Arguments$getReadablePathnames(reads2, absolutePath=TRUE)
    assertNoDuplicated(reads2)
  }

  # Argument 'indexPrefix':
  indexPrefix <- Arguments$getCharacter(indexPrefix, length=c(1L,1L))
  indexPath <- dirname(indexPrefix)
  indexName <- basename(indexPrefix)
  indexPath <- Arguments$getReadablePath(indexPath, absolutePath=TRUE)

  # Argument 'pathnameSAM':
  pathnameSAM <- Arguments$getWritablePathname(pathnameSAM, absolutePath=TRUE)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
     pushState(verbose)
     on.exit(popState(verbose), add=TRUE)
  }


  verbose && enter(verbose, "Running bowtie2()")

  verbose && cat(verbose, "R1 FASTQ files:")
  verbose && print(verbose, sapply(reads1, FUN=getRelativePath))
  if (isPaired) {
     verbose && cat(verbose, "R2 FASTQ files:")
     verbose && print(verbose, sapply(reads2, FUN=getRelativePath))
  }
  verbose && cat(verbose, "Paired alignment: ", isPaired)
  verbose && cat(verbose, "Reference index path: ", indexPath)
  verbose && cat(verbose, "Reference index prefix: ", indexPrefix)
  outPath <- dirname(pathnameSAM)
  verbose && cat(verbose, "Output directory: ", outPath)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Handle gzip'ed FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  pathnameFQ <- c(reads1, reads2)
  assertNoDuplicated(pathnameFQ)
  isGzipped <- any(sapply(pathnameFQ, FUN=isGzipped))
  if (isGzipped) {
    verbose && enter(verbose, "Detected gzipped FASTQ files")

    if (is.na(gzAllowed)) {
      gzAllowed <- queryBowtie2("support:fastq.gz")
    }

    if (!gzAllowed) {
      verbose && enter(verbose, "Temporarily decompressing gzipped FASTQ files")

      decompress <- getOption(aromaSettings, "devel/fastq.gz/decompress", TRUE)
      if (!decompress) {
        why <- attr(gzAllowed, "why")
        throw(sprintf("Cannot align reads in '%s': %s", pathnameFQ, why))
      }

      # If not, temporarily decompress (=remove when done)
      pathnameFQ <- sapply(pathnameFQ, FUN=gunzip, temporary=TRUE, remove=FALSE)
      pathnameFQtmpA <- pathnameFQ
      onExit({
        # Make sure to remove temporary file
        lapply(pathnameFQtmpA, FUN=function(pathname) {
          if (isFile(pathname)) file.remove(pathname)
        })
      })

      # Sanity check
      isGzipped <- any(sapply(pathnameFQ, FUN=isGzipped))
      .stop_if_not(!isGzipped)

      # Reassign 'reads1' and 'reads2'.
      idxs <- seq_along(reads1)
      reads1 <- pathnameFQ[ idxs]
      if (isPaired) reads2 <- pathnameFQ[-idxs]
      idxs <- NULL # Not needed anymore

      verbose && exit(verbose)
    }

    verbose && exit(verbose)
  } # if (isGzipped)
  pathnameFQ <- NULL # Not needed anymore



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate output atomically by writing to a temporary directory
  # that is renamed upon completion.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  outPathOrg <- outPath
  filenameSAM <- basename(pathnameSAM)
  filenameSAM <- Arguments$getTuxedoOption(filenameSAM)
  dirT <- tools::file_path_sans_ext(filenameSAM)
  outPath <- file.path(outPath, dirT)
  outPath <- sprintf("%s.TMP", outPath)
  outPath <- Arguments$getWritablePath(outPath, mustNotExist=TRUE)
  verbose && cat(verbose, "Temporary output directory: ", outPath)
  # At the end, assume failure, unless successful.
  outPathFinal <- sprintf("%s.ERROR", outPath)
  onExit({
    if (outPathFinal == outPathOrg) {
      # Move all files to destination directory
      files <- list.files(path=outPath, full.names=FALSE, recursive=FALSE)
      for (file in files) {
        src <- file.path(outPath, file)
        dest <- file.path(outPathFinal, file)
        file.rename(src, dest)
      }
      removeDirectory(outPath, recursive=FALSE, mustExist=FALSE)
    } else {
      # Tag working directory with an 'ERROR' tag
      file.rename(outPath, outPathFinal)
    }
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
  # Workaround the fact that bowtie2 binary does not support commas
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
  indexPath <- createLink(link=link, target=indexPath)
  onExit({ removeDirectory(indexPath) })
  link <- NULL # Not needed anymore
  indexPrefix <- file.path(indexPath, basename(indexPrefix))
  indexPrefix <- Arguments$getTuxedoOption(indexPrefix)

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

  # When reaching this point, 'bowtie' will not be able to tell whether
  # we are using the original files or the temporary workaround ones.
  # This is also reflected in the variable names, such that the below
  # code would look the same regardless whether we use a workaround or
  # not.  This means that if bowtie one day will support commas in
  # directories and filenames, we can just delete this whole section
  # and everything will work out of the box. /HB 2014-08-08


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call Bowtie2 executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && exit(verbose, "Calling systemBowtie2()")

  indexPrefix <- Arguments$getTuxedoOption(indexPrefix)
  args <- list("-x"=shQuote(indexPrefix))
  if (isPaired) {
    reads1 <- Arguments$getTuxedoOption(reads1)
    reads2 <- Arguments$getTuxedoOption(reads2)
    args[["-1"]] <- shQuote(paste(reads1, collapse=","))
    args[["-2"]] <- shQuote(paste(reads2, collapse=","))
  } else {
    reads1 <- Arguments$getTuxedoOption(reads1)
    args[["-U"]] <- shQuote(paste(reads1, collapse=","))
  }
  args[["-S"]] <- shQuote(filenameSAM)
  args <- c(args, list(...))

  verbose && cat(verbose, "Arguments:")
  verbose && print(verbose, args)

  res <- systemBowtie2(args=args, verbose=verbose)
  status <- attr(res, "status"); if (is.null(status)) status <- 0L
  verbose && cat(verbose, "Results:")
  verbose && str(verbose, res)
  verbose && cat(verbose, "Status:")
  verbose && str(verbose, status)

  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # In case bowtie2 generates empty SAM files
  # /HB 2012-10-01 (still to be observed)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (isFile(filenameSAM)) {
    if (file.info(filenameSAM)$size == 0L) {
      verbose && cat(verbose, "Removing empty SAM file falsely created by Bowtie2: ", filenameSAM)
      file.remove(filenameSAM)
      status <- -1L
    }
  }

  # Success?
  if (status == 0L && isFile(filenameSAM)) {
    outPathFinal <- outPathOrg
  }

  verbose && exit(verbose)

  res
} # bowtie2()
