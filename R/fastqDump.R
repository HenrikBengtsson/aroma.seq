###########################################################################/**
# @RdocDefault fastqDump
#
# @title "Calls the fastq-dump executable"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments passed to @see "systemFastqDump".}
#   \item{bSplit3}{If TRUE, use '--split-3' option to produce two fastq files for paired-end data}
#   \item{bGzip}{If TRUE, gzip output fastq files}
#   \item{outPath}{Directory where output fastq files will reside at completion}
#   \item{pathnames}{Zero or more input pathnames.}

#   \item{verbose}{See @see "R.utils::Verbose".}

# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("fastqDump", "default", function(...,
                                             bSplit3=TRUE,
                                             bGzip=TRUE,
                                             outPath="fastqData",
                                             pathnames=character(0L), verbose=FALSE) {
  
  # Required statement for onExit to be evaluated
  on.exit(eval(onExit()))
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  # Argument 'bSplit3':
  bSplit3 <- Arguments$getLogical(bSplit3)
  
  # Argument 'bGzip':
  bGzip <- Arguments$getLogical(bGzip)
  
  # Argument 'outPath'
  outPath <- Arguments$getWritablePath(outPath)
  
  # Argument 'pathnames':
  pathnames <- lapply(pathnames, FUN=Arguments$getReadablePathname)
  pathnames <- unlist(pathnames, use.names=TRUE)
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Running FastqDump")
  verbose && cat(verbose, "Input pathnames:")
  verbose && print(verbose, pathnames)

  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Generate output atomically by writing to a temporary directory
  # that is renamed upon completion.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  outPathTmp <- sprintf("%s.TMP", outPath)
  outPathTmp <- Arguments$getWritablePath(outPathTmp, mustNotExist=TRUE)
  verbose && cat(verbose, "Temporary output directory: ", outPathTmp)
  # At the end, assume failure, unless successful (and outPathFinal is changed to outPath)
  outPathFinal <- sprintf("%s.ERROR", outPathTmp)
  onExit({
    file.rename(outPathTmp, outPathFinal)
  })
  
    
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # While running, let the output directory be the working directory.
  # Make sure the working directory is restored when exiting.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  opwd <- setwd(outPathTmp)
  verbose && cat(verbose, "Original working directory: ", opwd)
  onExit({
    verbose && cat(verbose, "Resetting working directory: ", opwd)
    setwd(opwd)
  })
  
  args <- list(...)
  if (length(pathnames) > 0L) {
    # Always quote pathnames (to handle spaces)
    pathnames <- sprintf('"%s"', pathnames)
    args <- c(list(pathnames), args)
  }
  if (bSplit3) {
    args <- c(args, "--split-3")    
  }

  verbose && cat(verbose, "Arguments:")
  verbose && str(verbose, args)
  args <- c(args,verbose=less(verbose, 10))
  res <- do.call(systemFastqDump, args=args)
  
  status <- attr(res, "status"); if (is.null(status)) status <- 0L
  verbose && cat(verbose, "Results:")
  verbose && str(verbose, res)
  verbose && cat(verbose, "Status:")
  verbose && str(verbose, status)
  verbose && exit(verbose)
  
  # Successful?
  if (status == 0L) {
    # If we get this far, assume it was all successful.
    fastqFiles <- findFiles(path=".", pattern="[.]fastq$", firstOnly=FALSE)
    if (bGzip) {
        fastqFiles <- sapply(fastqFiles, function(f) gzip(f), simplify=TRUE)
    }
    # Allow the temporary output path to be renamed to the
    # intended output path instead of the "error" one.
    outPathFinal <- outPath
  }
  
  verbose && exit(verbose)
  res

}) # fastqDump()
