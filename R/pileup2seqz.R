pileup2seqz <- function(pus, gc, sampleName, dataset, tags="seqz", organism, pathD=file.path("seqzData", fullname(dataset, tags), organism), ..., force=FALSE, verbose=FALSE) {
  use("sequenza")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pus':
  stopifnot(is.list(pus))
  stopifnot(length(pus) == 2)

  # Argument 'gc':
  gc <- Arguments$getInstanceOf(gc, "GcBaseFile")
  pathnameGC <- Arguments$getReadablePathname(getPathname(gc))

  ## Argument 'sampleName':
  sampleName <- Arguments$getCharacter(sampleName)

  ## Argument 'dataset':
  dataset <- Arguments$getCharacter(dataset)

  ## Argument 'tags':
  tags <- Arguments$getTags(tags)

  ## Argument 'organism':
  organism <- Arguments$getCharacter(organism)

  # Argument 'pathD':
  pathD <- Arguments$getWritablePath(pathD)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "pileup2seqz")
  verbose && cat(verbose, "Normal:")
  ns <- as.list(pus[[1]])
  names(ns) <- gsub(".*,chr=", "", names(ns))
  verbose && print(verbose, MPileupFileSet(ns))

  verbose && cat(verbose, "Tumor:")
  ts <- as.list(pus[[2]])
  names(ts) <- gsub(".*,chr=", "", names(ts))
  verbose && print(verbose, MPileupFileSet(ts))

  verbose && print(verbose, gc)

  verbose && cat(verbose, "Sample name: ", sampleName)
  verbose && cat(verbose, "Output path: ", pathD)

  cmd <- system.file("exec", "sequenza-utils.py", package="sequenza", mustWork=TRUE)
  verbose && cat(verbose, "Sequenza executable: ", cmd)
  stopifnot(isCapableOf(aroma.seq, "python"))
  
  chromosomes <- names(ns)
  chromosomes <- gsub(".*,chr=", "", chromosomes)
  nchrs <- length(chromosomes)
  verbose && printf(verbose, "Chromosomes: [%d] %s\n", nchrs, hpaste(chromosomes))
  stopifnot(!is.null(chromosomes))

  plist <- listenv()
  for (chr in chromosomes) {
    verbose && enterf(verbose, "Chromosome '%s'", chr)

    chrTag <- sprintf("chr=%s", chr)
    chrTag <- sprintf("chr=chr%s", chr)
    filenameD <- sprintf("%s,%s.seqz", sampleName, chrTag)
    pathnameD <- file.path(pathD, filenameD)
    verbose && cat(verbose, "Output pathname: ", pathnameD)
    pathnameDz <- sprintf("%s.gz", pathnameD)

    ## Future label 
    label <- sprintf("chr_%s-%s", chr, sampleName)
    
    if (!force && !isFile(pathnameD) && !isFile(pathnameDz)) {
      verbose && enter(verbose, "Running sequenza-utils.py pileup2seqz")
      n <- ns[[chr]]
      t <- ts[[chr]]
      stopifnot(!is.null(n), !is.null(t))
      if (inherits(n, "MPileupFile")) n <- getPathname(n)
      if (inherits(t, "MPileupFile")) t <- getPathname(t)
      n <- Arguments$getReadablePathname(n)
      t <- Arguments$getReadablePathname(t)
      args <- c("pileup2seqz", "-gc", pathnameGC, "-n", n, "-t", t)
      verbose && printf(verbose, "Call: %s %s\n", cmd, paste(args, collapse=" "))
      mstr(args)
      plist[[chr]] %<-% {
        pathnameT <- pushTemporaryFile(pathnameD)
        verbose && printf(verbose, "Call: %s %s\n", cmd, paste(args, collapse=" "))
        res <- system2(cmd, args=args, stdout=pathnameT)
        verbose && printf(verbose, "Return code: %s\n", res)
        if (res != 0) {
	  throw(sprintf("Error (return code %s) while running system command: %s %s", res, cmd, paste(args, collapse=" ")))
	} else if (file.size(pathnameT) == 0) {
	  throw(sprintf("Result file of %s is empty (0 bytes): %s", args[1], pathnameT))
        }
        popTemporaryFile(pathnameT)
        SeqzFile(gzip(pathnameD))
      } %label% label
      verbose && exit(verbose)
   } else {
      if (isFile(pathnameDz)) {
        verbose && cat(verbose, "Already processed. Skipping.")
        plist[[chr]] <- SeqzFile(pathnameDz)
      } else {
        verbose && cat(verbose, "Already processed, but will compress.")
        plist[[chr]] %<-% SeqzFile(gzip(pathnameD)) %label% label
      }
    }

    verbose && exit(verbose)
  } # for (chr ...)

  verbose && exit(verbose)

  ## Resolve and collect values as they are ready
  plist <- resolve(plist, value = TRUE)

  ## Return as a gzip'ed SeqzFileSet
  plist <- as.list(plist)
  plist <- SeqzFileSet(plist)
  
  plist
} # pileup2seqz()
