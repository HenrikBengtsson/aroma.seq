###########################################################################/**
# @RdocClass Bowtie2Alignment
#
# @title "The Bowtie2Alignment class"
#
# \description{
#  @classhierarchy
#
#  ...
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "AbstractAlignment".}
#  \item{groupBy}{A @character string or an explicit named @list,
#   specifying which input files should be processed together.}
#  \item{indexSet}{An @see "Bowtie2IndexSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Supported operating systems}{
#   This method is available on Linux, macOS, and Windows [1].
# }
#
# \author{Henrik Bengtsson and Pierre Neuvial}
#
# \references{
#  [1] Bowtie2, John Hopkins University, 2013.
#      \url{http://bowtie-bio.sourceforge.net/bowtie2/}
# }
#*/###########################################################################
setConstructorS3("Bowtie2Alignment", function(..., groupBy=NULL, indexSet=NULL) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'groupBy':
  if (is.null(groupBy)) {
  } else if (is.character(groupBy)) {
    groupBy <- match.arg(groupBy, choices=c("name"))
  } else if (is.list(groupBy)) {
    # Validated below
  } else {
    throw("Invalid argument 'groupBy': ", mode(groupBy))
  }

  # Argument 'indexSet':
  if (!is.null(indexSet)) {
    indexSet <- Arguments$getInstanceOf(indexSet, "Bowtie2IndexSet")
  }

  # Arguments '...':
  args <- list(...)

  this <- extend(AbstractAlignment(..., indexSet=indexSet, groupBy=groupBy), c("Bowtie2Alignment", uses("FileGroupsInterface")))

  # Argument 'groupBy':
  if (is.list(groupBy)) {
    validateGroups(this, groups=groupBy)
  }

  this
})


setMethodS3("getSampleNames", "Bowtie2Alignment", function(this, ...) {
  getGroupNames(this, ...)
}, protected=TRUE)


###########################################################################/**
# @RdocMethod process
#
# @title "Runs the aligner"
#
# \description{
#   @get "title" on all input files.
#   The generated BAM files are all sorted and indexed.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
#  \item{skip}{If @TRUE, already processed files are skipped.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "BamDataSet".
# }
#
# @author
#*/###########################################################################
setMethodS3("process", "Bowtie2Alignment", function(this, ..., skip=TRUE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  asBowtie2Parameters <- function(rg, ...) {
    if (isEmpty(rg)) {
      return(NULL)
    }

    # Validate
##    if (!hasID(rg)) {
##      throw("Bowtie requires that the SAM read group has an ID.")
##    }

    rgId <- asSamList(rg)$ID
    if (is.null(rgId)) rgId <- 1L

    rgArgs <- asString(rg, fmtstr="%s:%s")
    rgArgs <- rgArgs[regexpr("^ID:", rgArgs) == -1L]

    # Don't forget to put within quotation marks
    rgArgs <- sprintf("\"%s\"", rgArgs)

    # Escape spaces [NOT ENOUGH, because bowtie2 calls 'perl bowtie2'
    # and in the latter step it is all lost. /HB 2014-08-11]
    ## rgArgs <- gsub(" ", "\\ ", rgArgs, fixed=TRUE)

    # Sanity check
    nok <- hasCommas(rgArgs)
    if (any(nok)) {
      throw("SAM Read Group options must not contain spaces (=not supported by Bowtie2): ", hpaste(sQuote(rgArgs[nok])))
    }

    rgArgs <- as.list(rgArgs)
    names(rgArgs) <- rep("--rg", times=length(rgArgs))

    rgArgs <- c(list("--rg-id"=rgId), rgArgs)

    rgArgs
  } # asBowtie2Parameters()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'force':
  force <- Arguments$getLogical(force)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Bowtie2 alignment")
  ds <- getInputDataSet(this)
  verbose && cat(verbose, "Input data set:")
  verbose && print(verbose, ds)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Get groups of items to be processed at the same time
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Grouping input data set")
  groups <- getGroups(this)
  verbose && printf(verbose, "Merging into %d groups: %s\n", length(groups), hpaste(names(groups)))
  verbose && str(verbose, head(groups))
  verbose && cat(verbose, "Number of items per groups:")
  ns <- sapply(groups, FUN=length)
  t <- table(ns)
  names(t) <- sprintf("n=%s", names(t))
  verbose && print(verbose, t)
  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify groups to be processed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (force) {
    todo <- seq_along(groups)
  } else {
    bams <- getOutputDataSet(this, onMissing="NA", verbose=less(verbose, 1))
    todo <- which(!sapply(bams, FUN=isFile))
  }
  verbose && cat(verbose, "Number of groups to process: ", length(todo))

  # Already done?
  if (!force && length(todo) == 0L) {
    verbose && cat(verbose, "Already processed.")
    verbose && print(verbose, bams)
    verbose && exit(verbose)
    return(bams)
  }

  isPaired <- isPaired(ds)
  verbose && cat(verbose, "Paired-end analysis: ", isPaired)

  outPath <- getPath(this)
  verbose && cat(verbose, "Output directory: ", outPath)

  # Additional alignment parameters
  params <- getParameters(this)
  verbose && cat(verbose, "Parameters:")
  verbose && str(verbose, params)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Retrieving/building Bowtie2 index set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Retrieving/building Bowtie2 index set")
  is <- getIndexSet(this)
  verbose && cat(verbose, "Aligning using index set:")
  verbose && print(verbose, is)
  indexPrefix <- getIndexPrefix(is)
  verbose && exit(verbose)

  rgSet <- this$.rgSet
  if (!is.null(rgSet)) {
    verbose && cat(verbose, "Assigning SAM read group:")
    verbose && print(verbose, rgSet)
    validate(rgSet)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # User arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- params
  # Drop already used parameters
  args$groupBy <- NULL
  verbose && cat(verbose, "User arguments:")
  verbose && str(verbose, args)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Number of files: ", length(ds))
  verbose && cat(verbose, "Number of groups: ", length(groups))

  isPaired <- isPaired(this)
  outPath <- getPath(this)

  res <- listenv()
  IDXS <- groups[todo]
  for (gg in seq_along(IDXS)) {
    idxs <- IDXS[[gg]]
    dfListR1 <- as.list(ds[idxs])
    sampleName <- names(IDXS)[gg]

    verbose && enter(verbose, sprintf("Bowtie2 alignment sample #%d (%s) of %d", gg, sampleName, length(IDXS)))

    # The BAM file to be generated
    fullname <- sampleName
    filename <- sprintf("%s.bam", fullname)
    pathnameBAM <- Arguments$getWritablePathname(filename, path=outPath)
    verbose && cat(verbose, "BAM pathname: ", pathnameBAM)

    done <- (skip && isFile(pathnameBAM))
    if (done) {
      verbose && cat(verbose, "Already aligned. Skipping")
      res[[gg]] <- pathnameBAM
      verbose && exit(verbose)
      next
    }

    res[[gg]] %<-% {
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (a) Generate SAM file
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      filename <- sprintf("%s.sam", fullname)
      pathnameSAM <- Arguments$getWritablePathname(filename, path=outPath)
      verbose && cat(verbose, "Temporary/intermediate SAM pathname: ", pathnameSAM)
      if (!isFile(pathnameSAM)) {
        reads1 <- sapply(dfListR1, FUN=getPathname)
        verbose && printf(verbose, "R1 FASTQ files: [%d] %s\n", length(reads1), hpaste(sQuote(reads1)))

        # Final sample-specific output directory
        args <- list(
          indexPrefix=indexPrefix,
          reads1=reads1,
          reads2=NULL,
          ...,
          pathnameSAM=pathnameSAM
        )

        if (isPaired) {
          dfListR2 <- lapply(dfListR1, FUN=getMateFile)
          reads2 <- sapply(dfListR2, FUN=getPathname)
          verbose && printf(verbose, "R2 FASTQ files: [%d] %s\n", length(reads2), hpaste(sQuote(reads2)))
          args$reads2 <- reads2
        }

        # Extract sample-specific read group
        rgII <- getSamReadGroup(dfListR1[[1L]])
        if (length(rgSet) > 0L) {
          rgII <- merge(rgSet, rgII)
        }
        verbose && cat(verbose, "Writing SAM Read Groups:")
        verbose && print(verbose, rgII)
        verbose && cat(verbose, "Bowtie2 parameters:")
        rgArgs <- asBowtie2Parameters(rgII)
        verbose && print(verbose, rgArgs)
        rgII <- NULL  # Not needed anymore

        args <- c(args, rgArgs)
        verbose && cat(verbose, "Arguments:")
        verbose && str(verbose, args)
        args$verbose <- less(verbose, 5)

        res <- do.call(bowtie2, args=args)
        verbose && cat(verbose, "System result code: ", res)
      } # if (!isFile(pathnameSAM))

      # Sanity check
      Arguments$getReadablePathname(pathnameSAM)

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (b) Generates a (sorted and indexed) BAM file from SAM file
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!isFile(pathnameBAM)) {
        sf <- SamDataFile(pathnameSAM)
        bf <- convertToBam(sf, verbose=less(verbose, 5))
        verbose && print(verbose, pathnameBAM)
      }
      # Sanity check
      Arguments$getReadablePathname(pathnameBAM)

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (c) Remove temporary/intermediate SAM file
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      file.remove(pathnameSAM)

      pathnameBAM
    } ## %<-%

    verbose && exit(verbose)
  } ## for (gg ...)

  ## Resolve futures
  res <- resolve(res)

  bams <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1))
  verbose && print(verbose, bams)

  ## Sanity checks
  is <- getIndexSet(this)
  for (ii in seq_along(bams)) {
    bam <- bams[[ii]]
    .stop_if_not(isFile(bam))
    isCompatibleWith(bam, is)
  }

  verbose && exit(verbose)

  bams
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# TO DROP
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("validateGroups", "Bowtie2Alignment", function(this, groups, ...) {
  # Input data set
  ds <- getInputDataSet(this)
  nbrOfFiles <- length(ds)

  # Sanity checks
  idxs <- unlist(groups, use.names=FALSE)
  idxs <- Arguments$getIndices(idxs, max=nbrOfFiles)
  if (length(idxs) < nbrOfFiles) {
    throw("One or more input FASTQ files is not part of any group.")
  } else if (length(idxs) > nbrOfFiles) {
    throw("One or more input FASTQ files is part of more than one group.")
  }

  if (is.null(names(groups))) {
    throw("The list of groups does not have names.")
  }

  invisible(groups)
}, protected=TRUE)
