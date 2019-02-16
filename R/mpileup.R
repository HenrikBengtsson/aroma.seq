#' Mpileup on one BAM file over one or more chromosomes
#'
#' @param bams A BamDataFileSet.
#' @param fa A FastaReferenceFile.
#' @param Q A numeric quality threshold.
#' @param chromosomes A character vector of chromosomes to process.
#' @param pathD The directory where to write the mpileup files.
#' @param ... Passed to \code{mpileup()} for BamDataFile.
#' @param force If TRUE, already processes samples and / or chromosomes are reprocessed.
#' @param verbose If TRUE, verbose output is printed.
#'
#' @return A list of lists of gzip'ed MpileupFile:s.
#'
#' @export
setMethodS3("mpileup", "BamDataSet", function(bams, fa, Q=20, chromosomes=getSeqNames(fa), dataset=getFullName(bams), tags="mpileup", organism=getOrganism(bams), pathD=file.path("seqzData", fullname(dataset, tags), organism), ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'bams':
  .stop_if_not(inherits(bams, "BamDataSet"))
  .stop_if_not(isCompatibleWith(bams, bams[[1]]))  ## Self consistency

  # Argument 'fa':
  .stop_if_not(inherits(fa, "FastaReferenceFile"))
  ## Assert compatible chromosome names and equal chromosome lengths
  seqLengths <- getSeqLengths(fa)
  for (ii in seq_along(bams)) {
    bam <- bams[[ii]]
    .stop_if_not(all(getSeqNames(bam) %in% names(seqLengths)))
    .stop_if_not(all(seqLengths[getSeqNames(bam)] == getSeqLengths(bam)))
  }
  
  # Argument 'pathD':
  pathD <- Arguments$getWritablePath(pathD)

  # Argument 'chromosomes':
  chromosomes <- Arguments$getCharacters(chromosomes)
  if (any(duplicated(chromosomes))) {
    throw("Argument 'chromosomes' contains duplicates: %s", hpaste(unique(chromosomes[duplicated(chromosomes)])))
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "mpileup()")
  verbose && print(verbose, bams)
  verbose && print(verbose, fa)
  verbose && cat(verbose, "Output path: ", pathD)

  nchrs <- length(chromosomes)
  verbose && printf(verbose, "Chromosomes: [%d] %s\n", nchrs, hpaste(chromosomes))

  ## Assert BAM use the same chromosome names
  .stop_if_not(all(chromosomes %in% getSeqNames(fa)))

  ## Assert BAM are indexed, i.e. they have BAI files
  .stop_if_not(all(vapply(bams, FUN = hasIndex, FUN.VALUE = FALSE)))

  ## Chromosome-by-sample matrix of MPileupFile:s
  ## (will be transposed before being returned)
  res <- list()
  dimnames <- list(chromosomes, getFullNames(bams))
  dim <- lengths(dimnames)
  length(res) <- prod(dim)
  dim(res) <- dim
  dimnames(res) <- dimnames
 
  ## Step 1: Grab already finished results
  if (!force) {
    verbose && enter(verbose, "Identify which BAM files and chromosomes are to be processed")
    fullnames <- getFullNames(bams)
    chrTags <- sprintf("chr=%s", gsub("^chr", "", chromosomes))
    filenamesDz <- sprintf("%s,%s.mpileup.gz", rep(fullnames, each=length(chromosomes)), chrTags)
    pathnamesDz <- file.path(pathD, filenamesDz)
    done <- which(file_test("-f", pathnamesDz))
    ndone <- length(done)
    verbose && printf(verbose, "Already processed: %d (%.1f%%) out of %d (%dx%d)\n",
           ndone, 100*ndone/length(pathnamesDz), length(pathnamesDz), length(bams), nchrs)
    if (ndone > 0) {
      verbose && enterf(verbose, "Setting up %d existing output files", ndone)
      res[done] <- lapply(pathnamesDz[done], FUN=MPileupFile)
      verbose && exit(verbose)
    }
    verbose && print(verbose, res)
    verbose && exit(verbose)
  } # if (!force)


  ## Step 2: Process incomplete samples (i.e. missing one or more chromosomes)
  verbose && enter(verbose, "Processing incomplete BAMs (i.e. missing one or more chromosomes)")

  work <- listenv()
  length(work) <- length(bams)
  for (ii in seq_along(bams)) {
    bam <- bams[[ii]]
    name <- getFullName(bam)
    verbose && enterf(verbose, "BAM #%d ('%s') of %d", ii, name, length(bams))

    ## Already processed?
    resII <- res[,ii]
    todo <- which(unlist(lapply(resII, FUN=is.null)))
    if (length(todo) == 0) {
      verbose && cat(verbose, "All chromosomes already processed. Skipping.")
      verbose && exit(verbose)
      next
    }
    
    verbose && print(verbose, bam)
    verbose && print(verbose, fa)
    chromosomesT <- chromosomes[todo]
    nchrsT <- length(chromosomesT)
    verbose && printf(verbose, "Chromosomes to process: [%d] %s\n", nchrsT, hpaste(chromosomesT))

    ## Assert BAM use the same chromosome names
    .stop_if_not(all(chromosomesT %in% getSeqNames(bam)))

    label <- sprintf("sample_%d", ii)
    work[[ii]] %<-% {
      print(future::sessionDetails())
      resT <- mpileup(bam, fa=fa, chromosomes=chromosomesT, Q=Q, pathD=pathD, ..., force=force, verbose=less(verbose, 1))
      verbose && print(verbose, resT)
      .stop_if_not(length(resT) == length(todo))
      resT
    } %label% label
    
    verbose && exit(verbose)
  } # for (ii ...)
  todo <- NULL
  verbose && print(verbose, work)
  verbose && exit(verbose)


  ## Step 3: Resolve and collect
  verbose && enter(verbose, "Collect and resolve all futures")
  work <- resolve(work, value=TRUE)
  verbose && print(verbose, work)
  .stop_if_not(length(work) == length(bams))
  verbose && exit(verbose)

  ## Step 4: Gather and rearrange
  verbose && enter(verbose, "Gather and rearrange")
  .stop_if_not(length(work) == ncol(res))
  for (ii in seq_along(work)) {
    resII <- res[,ii]
    todo <- which(unlist(lapply(resII, FUN=is.null)))
    if (length(todo) == 0) next
    workII <- work[[ii]]
    .stop_if_not(length(workII) == length(todo))
    res[todo, ii] <- workII
  }
  work <- todo <- NULL ## Not needed anymore

  ## Transpose
  res <- t(res)
  verbose && str(verbose, res)
  verbose && exit(verbose)

  verbose && exit(verbose)

  res
}) # mpileup() for BamDataSet


#' @param bam A BamDataFile.
#'
#' @return A list of \code{length(chromosomes)} gzip'ed MPileupFile:s.
#'
#' @export
setMethodS3("mpileup", "BamDataFile", function(bam, fa, Q=20, chromosomes=getSeqNames(fa), pathD, ..., force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'bam':
  .stop_if_not(inherits(bam, "BamDataFile"))

  # Argument 'fa':
  .stop_if_not(inherits(fa, "FastaReferenceFile"))
  ## Assert compatible chromosome names and equal chromosome lengths
  seqLengths <- getSeqLengths(fa)
  .stop_if_not(all(getSeqNames(bam) %in% names(seqLengths)))
  .stop_if_not(all(seqLengths[getSeqNames(bam)] == getSeqLengths(bam)))
  
  # Argument 'pathD':
  pathD <- Arguments$getWritablePath(pathD)

  # Argument 'chromosomes':
  chromosomes <- Arguments$getCharacters(chromosomes)
  if (any(duplicated(chromosomes))) {
    throw("Argument 'chromosomes' contains duplicates: %s", hpaste(unique(chromosomes[duplicated(chromosomes)])))
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "mpileup()")
  verbose && print(verbose, bam)
  verbose && print(verbose, fa)
  verbose && cat(verbose, "Output path: ", pathD)

  nchrs <- length(chromosomes)
  verbose && printf(verbose, "Chromosomes: [%d] %s\n", nchrs, hpaste(chromosomes))

  ## Assert BAM use the same chromosome names
  .stop_if_not(all(chromosomes %in% getSeqNames(fa)))
  .stop_if_not(all(chromosomes %in% getSeqNames(bam)))
  
  name <- getFullName(bam)

  res <- listenv()
  names(res) <- names(chromosomes)
  
  for (kk in seq_along(chromosomes)) {
    chr <- chromosomes[kk]
    verbose && enterf(verbose, "Chromosome '%s'", chr)

    chrTag <- sprintf("chr=%s", gsub("^chr", "", chr))
    filenameD <- sprintf("%s,%s.mpileup", name, chrTag)
    pathnameD <- file.path(pathD, filenameD)
    pathnameDz <- sprintf("%s.gz", pathnameD)
    verbose && cat(verbose, "Output pathname: ", pathnameDz)

    ## Future label
    label <- sprintf("chr_%s", chr)

    if (!force && !isFile(pathnameD) && !isFile(pathnameDz)) {
      verbose && enter(verbose, "Running samtools mpileup")
      res[[kk]] %<-% {
        print(future::sessionDetails())
        pathnameT <- pushTemporaryFile(pathnameD)
	ans <- samtoolsMpileup(refFile=getPathname(fa), bamFile=getPathname(bam), pathnameD=pathnameT, r=chr, Q=Q, verbose=verbose)
	if (ans == 0L) popTemporaryFile(pathnameT)
        mp <- MPileupFile(pathnameD)
	gzip(mp)
	mp
     } %label% label
      verbose && exit(verbose)
    } else {
      if (isFile(pathnameDz)) {
        verbose && cat(verbose, "Already processed. Skipping.")
        res[[kk]] <- MPileupFile(pathnameDz)
      } else {
        verbose && cat(verbose, "Already processed, but will compress.")
        res[[kk]] %<-% {
          gzip(pathnameD)
	  MPileupFile(pathnameDz)
	} %label% label
      }
    }

    verbose && exit(verbose)
  } # for (kk ...)

  verbose && enter(verbose, "Collect and resolve all futures")
  res <- resolve(res, value=TRUE)
  verbose && print(verbose, res)
  verbose && exit(verbose)

  ## Coerce to a list
  res <- as.list(res)
  verbose && str(verbose, res)

  verbose && exit(verbose)

  res
}) # mpileup() for BamDataFile
