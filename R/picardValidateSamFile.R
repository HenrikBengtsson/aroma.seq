setMethodS3("picardValidateSamFile", "default", function(pathname, ..., onWarning=c("warning", "error", "ignore"), onError=c("error", "warning", "ignore"), verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local function
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  reportOn <- function(..., onAction=c("warning", "error", "ignore")) {
    onAction <- match.arg(onAction)
    msg <- paste(..., collapse="")
    if (onAction == "error") {
      throw(msg)
    } else if (onAction == "warning") {
      warning(msg)
    }
  } # reportOn()

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  pathname <- Arguments$getReadablePathname(pathname)

  # Argument 'onWarning' and 'onError':
  onWarning <- match.arg(onWarning)
  onError <- match.arg(onError)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Validating SAM/BAM file using Picard")
  verbose && cat(verbose, "Pathname: ", pathname)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Run 'picard'
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Running 'picard ValidateSamFile'")
  path <- tempfile()
  filenameR <- sprintf("%s.log", basename(pathname))
  pathnameR <- Arguments$getWritablePathname(filenameR, path=path)
  on.exit({
    file.remove(pathnameR)
  }, add=TRUE)
  verbose && cat(verbose, "Temporary log file: ", pathnameR)

  res <- systemPicard("ValidateSamFile", INPUT=pathname, ..., stdout=pathnameR, stderr=pathnameR, verbose=less(verbose, 10))
  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Parse results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Parsing results")
  bfr <- readLines(pathnameR)
  for (type in c("WARNING", "ERROR")) {
    pattern <- sprintf("^%s: ", type)
    values <- grep(pattern, bfr, value=TRUE)
    if (length(values) == 0L) next

    verbose && print(verbose, values)
    msgs <- gsub(pattern, "", values)
    msgs <- paste(msgs, collapse="\n")
    if (type == "WARNING") {
      reportOn(msgs, onAction=onWarning)
    } else if (type == "ERROR") {
      reportOn(msgs, onAction=onError)
    }
  } # for (type ...)
  verbose && exit(verbose)

  verbose && exit(verbose)

  invisible(bfr)
}, protected=TRUE)
