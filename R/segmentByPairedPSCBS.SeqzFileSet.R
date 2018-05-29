setMethodS3("readTotalCNsAndBAFs", "SeqzFile", function(this, ploidy=2, ..., verbose=FALSE) {
  use("readr")
  
  verbose <- Arguments$getVerbose(verbose)
  pathname <- getPathname(this)

  col_types <- readr::cols_only(
    chromosome      = readr::col_character(),
    position        = readr::col_integer(),
    depth.normal    = readr::col_integer(),
    depth.tumor     = readr::col_integer(),
    Af              = readr::col_double(),  ## FIXME: Why not just void this one? /HB 2017-06-26
    Bf              = readr::col_double(),
    zygosity.normal = readr::col_character()
  )
  data <- readr::read_tsv(pathname, col_types=col_types, progress=isVisible(verbose))

  depth <- data$depth.tumor + data$depth.normal
  tcn <- ploidy * data$depth.tumor/data$depth.normal

  data <- data.frame(
    chromosome=data$chromosome, x=data$position,
    depth=depth, tcn=tcn, baf=data$Bf, isHet=(data$zygosity.normal == "het"), 
    stringsAsFactors=FALSE
  )

  data
})


setMethodS3("readStats", "SeqzFile", function(this, ..., force=FALSE, verbose=FALSE) {
  use("readr")
  
  verbose <- Arguments$getVerbose(verbose)

  fullname <- getFullName(this)
  filename <- sprintf("%s.seqz.stats.rds", fullname)
  pathnameS <- Arguments$getReadablePathname(filename, path=getPath(this), mustExist=FALSE)
  if (!force && isFile(pathnameS)) {
    stats <- loadObject(pathnameS)
    return(stats)
  }

  pathnameS <- Arguments$getWritablePathname(pathnameS, mustNotExist=FALSE)
  
  pathname <- getPathname(this)

  col_types <- readr::cols_only(
    depth.normal = readr::col_integer(),
    depth.tumor  = readr::col_integer()
  )
  data <- readr::read_tsv(pathname, col_types=col_types, progress=isVisible(verbose))

  stats <- apply(data, MARGIN=2L, FUN=median, na.rm=TRUE)
  stats <- c(stats, count=nrow(data))

  rm(list="data") ## Not needed anymore

  saveObject(stats, file=pathnameS)

  stats
})


setMethodS3("readStats", "SeqzFileSet", function(this, ..., verbose=FALSE) {
  verbose <- Arguments$getVerbose(verbose)

  verbose && enterf(verbose, "Collection read statistics for %d Sequenza files", length(this))

  stats <- list()
  for (ii in seq_along(this)) {
    df <- this[[ii]]
    verbose && enterf(verbose, "File #%d ('%s') of %d", ii, getFullName(df), length(this))
    stats[[ii]] <- readStats(df, ..., verbose=less(verbose, 10))
    verbose && exit(verbose)
  }
  names(stats) <- names(this)

  verbose && exit(verbose)

  stats <- Reduce(rbind, stats)
  rownames(stats) <- names(this)

  stats
})


setMethodS3("readDataFrameForPSCBS", "SeqzFile", function(this, ..., verbose=FALSE) {
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enterf(verbose, "Loading %s data", class(this)[1])
  
  data <- readTotalCNsAndBAFs(this, ..., verbose=less(verbose))
  verbose && str(verbose, data)

  ## Convert chromosome IDs to indices
  chr <- data$chromosome
  chr <- gsub("chr", "", chr, fixed=TRUE)
  chr[chr == "X"] <- 23
  chr[chr == "Y"] <- 24
  chr[chr == "M"] <- 25
  chr <- as.integer(chr)
  data$chromosome <- chr
  
  data$rho <- 2 * abs(data$baf - 1/2)
  data$rho[!data$isHet] <- NA_real_
  verbose && str(verbose, data)

  verbose && exit(verbose)

  data
})

setMethodS3("binForPSCBS", "data.frame", function(data, binSize, ..., verbose=FALSE) {
  # Argument 'binSize':
  binSize <- Arguments$getNumeric(binSize, range=c(1,Inf))
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Binning total and DH")
  chrs <- sort(unique(data$chromosome))
  verbose && cat(verbose, "Chromosomes: ", hpaste(chrs))

  dataS <- list()
  for (kk in seq_along(chrs)) {
    chr <- chrs[kk]
    verbose && enter(verbose, sprintf("Chromosome #%d ('%s') of %d", kk, chr, length(chrs)))
    
    dataT <- data[data$chromosome == chr,]
    
    xr <- range(dataT$x, na.rm=TRUE) / 1e6
    verbose && printf(verbose, "Positional range: [%g,%g] Mb\n", xr[1], xr[2])
    bx <- seq(from=floor(xr[1]), to=ceiling(xr[2]), by=binSize/1e6)*1e6
    verbose && print(verbose, "Bins boundaries:")
    verbose && str(verbose, bx)
    verbose && print(verbose, "Bins centers:")
    x <- (bx[-length(bx)] + bx[-1]) / 2
    verbose && str(verbose, x)

    depth <- binMeans(y=dataT$depth, w=dataT$depth, x=dataT$x, bx=bx, na.rm=TRUE)
    tcn <- binMeans(y=dataT$tcn, w=dataT$depth, x=dataT$x, bx=bx, na.rm=TRUE)
    rho <- binMeans(y=dataT$rho, w=dataT$depth, x=dataT$x, bx=bx, na.rm=TRUE)
    keep <- (depth > 0) & (is.finite(tcn) | is.finite(rho))
    x <- x[keep]
    tcn <- tcn[keep]
    rho <- rho[keep]
    depth <- depth[keep]
    dataT <- data.frame(chromosome=chr, x=x, CT=tcn, rho=rho, depth=depth, stringsAsFactors=FALSE)
    verbose && str(verbose, dataT)
    x <- tcn <- rho <- NULL

    dataS[[kk]] <- dataT
    
    verbose && exit(verbose)
  } ## for (kk ...)

  dataS <- Reduce(rbind, dataS)
  verbose && exit(verbose)
  
  dataS
}) ## binForPSCBS()



setMethodS3("segmentByPairedPSCBS", "SeqzFileSet", function(seqz, ..., binSize, minLength=2e6, dataset=getFullName(seqz), tags=sprintf("%gkb", binSize/1e3), organism=getSubdirs(seqz), outPath=file.path(rootPath, paste(c(dataset, tags), collapse=","), organism), rootPath="pscbsData", verbose=FALSE) {
  use("PSCBS (>= 0.61.0)")     ## segmentByPairedPSCBS()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'binSize':
  binSize <- Arguments$getNumeric(binSize, range=c(1,Inf))

  # Argument 'minLength':
  minLength <- Arguments$getNumeric(minLength, range=c(0,Inf))

  # Argument 'outPath':
  outPath <- Arguments$getWritablePath(outPath)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enterf(verbose, "Paired PSCBS on %s", class(seqz)[1])
  verbose && print(verbose, seqz)
  names <- unique(getNames(seqz))
  verbose && printf(verbose, "Samples: [%d] %s", length(names), hpaste(names))
  chrs <- unique(grep("^chr", sapply(seqz, FUN=getTags), value=TRUE))
  chrs <- gsub("chr(|=)", "", chrs)
  verbose && printf(verbose, "Chromosomes: [%d] %s\n", length(chrs), hpaste(chrs))
  verbose && printf(verbose, "Bin size: %g (%g kb)\n", binSize, binSize/1e3)

  verbose && cat(verbose, "Output path: ", outPath)

  fits <- listenv()
  for (ii in seq_along(seqz)) {
    seq <- seqz[[ii]]
    sample <- getFullName(seq)
    sample <- gsub(",chr(|=).*", "", sample)
    chrT <- gsub("chr(|=)", "", grep("^chr", getTags(seq), value=TRUE))
    chr <- as.integer(chrT)
    chr[chrT == "X"] <- 23L
    chr[chrT == "Y"] <- 24L
    chr[chrT == "M"] <- 25L
    verbose && enterf(verbose, "Chromosome #%d ('%s') of %d", ii, chr, length(seqz))
    verbose && printf(verbose, "Sample '%s'; Chromosome %d\n", sample, chr)

    chrTag <- sprintf("chr=%d", chr)

    ## Future label
    label <- sprintf("%s-%s", chrTag, sample)
    
    filenameD <- sprintf("%s,%s,PairedPSCBS.xdr", sample, chrTag)
    pathnameD <- file.path(outPath, filenameD)
    if (!isFile(pathnameD)) {
      fits[[chrTag]] %<-% {
        verbose <- Arguments$getVerbose(-10, timestamp=TRUE)

        verbose && enter(verbose, "Reading seqz data")
        verbose && print(verbose, seq)
        # Load data
        data <- readDataFrameForPSCBS(seq, verbose=less(verbose))
        verbose && str(verbose, data)
        nbrOfLoci <- nrow(data)
        nbrOfHets <- sum(data$isHet)
        verbose && exit(verbose)

        verbose && enter(verbose, "Binning total and DH")
	data <- binForPSCBS(data, binSize=binSize, verbose=less(verbose))
        verbose && str(verbose, data)
        verbose && exit(verbose)

        gc()

        verbose && enter(verbose, "Segmenting PSCN")
        if (minLength > 0) {
          gaps <- PSCBS::findLargeGaps(data, minLength=minLength)
          knownSegments <- PSCBS::gapsToSegments(gaps)
        } else {
          knownSegments <- NULL
        }

        fit <- PSCBS::segmentByPairedPSCBS(data, knownSegments=knownSegments, verbose=verbose)
        ## Ideally sampleName(fit) <- sample but 'future' complains!
        fit <- PSCBS::setSampleName(fit, sample)
        saveObject(fit, file=pathnameD)
        data <- NULL
        verbose && exit(verbose)

        verbose && print(verbose, fit)
        fit
      } %label% label
    } else {
      verbose && cat(verbose, "Already processed. Skipping. Loading.")
      fit <- loadObject(pathnameD)
      fits[[chrTag]] <- fit
    }
    verbose && exit(verbose)
  } # for (ii ...)

  ## Resolve futures and gather their values
  fits <- resolve(fits, value = TRUE)

  ## Coerce to a list
  fits <- as.list(fits)
  
  verbose && print(verbose, fits)
  verbose && exit(verbose)

  fits
}) # segmentByPairedPSCBS()
