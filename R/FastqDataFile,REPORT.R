setMethodS3("report", "FastqDataFile", function(this, dataSet, ..., flavor="qrqc", type="md", outPath=".", verbose=FALSE) {
  use("R.rsp")
  rfile <- R.rsp::rfile


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet)

  # Argument 'flavor':
  flavor <- match.arg(flavor)

  # Argument 'outPath':
  outPath <- Arguments$getWritablePath(outPath)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Generating report")
  verbose && print(verbose, this)
  verbose && cat(verbose, "Data set: ", dataSet)
  verbose && cat(verbose, "Flavor: ", flavor)
  verbose && cat(verbose, "Output type: ", type)

  # Locate report template
  pathname <- findRspReportTemplate(this, tags=flavor, type=type)
  if (length(pathname) == 0L) {
    throw("Failed to locate report template for this flavor and type: ", paste(c(sQuote(flavor), sQuote(type)), collapse=" & "))
  }
  verbose && cat(verbose, "Report template (source): ", pathname)

  # Output path
  outPath <- file.path(outPath, getFullName(this))
  outPath <- Arguments$getWritablePath(outPath)
  verbose && cat(verbose, "Report path: ", outPath)

  # Rename template
  filename <- basename(pathname)
  filenameT <- gsub("^[^,]*", getFullName(this), filename)
  pathnameT <- file.path(outPath, filenameT)
  copyFile(pathname, pathnameT, copy.mode=FALSE, overwrite=TRUE)
  verbose && cat(verbose, "Report template (renamed): ", pathnameT)

  # Generate PDF report
  rspArgs <- list(this, dataSet=dataSet)
  pathnameR <- rfile(pathnameT, ..., workdir=outPath, verbose=less(verbose, 5))

  verbose && cat(verbose, "Generated report: ", pathnameR)

  verbose && exit(verbose)

  invisible(pathnameR)
}) # report()


setMethodS3("report", "FastqDataSet", function(this, dataSet=getFullName(this), ..., flavor="qrqc", type="md", outPath=file.path("reports", fullname(dataSet, tags=flavor)), verbose=FALSE) {
  use("R.rsp")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet)

  # Argument 'flavor':
  flavor <- match.arg(flavor)

  # Argument 'outPath':
  outPath <- Arguments$getWritablePath(outPath)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Generating reports")
  verbose && print(verbose, this)
  verbose && cat(verbose, "Flavor: ", flavor)
  verbose && cat(verbose, "Output type: ", type)

  # Locate report template
  pathname <- findRspReportTemplate(this, tags=flavor, type=type)
  if (length(pathname) == 0L) {
    throw("Failed to locate report template for this flavor and type: ", paste(c(sQuote(flavor), sQuote(type)), collapse=" & "))
  }
  verbose && cat(verbose, "Report template (source): ", pathname)
  verbose && cat(verbose, "Output path: ", outPath)

  ## FIXME: report(df, ...) does not work with multiprocess futures
  ##        /HB 2016-03-12
  oplan <- future::plan("eager")
  on.exit(future::plan(oplan), add=TRUE)

  res <- listenv()
  for (ii in seq_along(this)) {
    df <- this[[ii]]
    name <- getFullName(df)
    verbose && enter(verbose, sprintf("FASTQ report #%d ('%s') of %d", ii, name, length(this)))

    res[[ii]] %<-% {
      report(df, dataSet=dataSet, flavor=flavor, outPath=outPath, ..., verbose=verbose)
    }

    verbose && exit(verbose)
  } ## for (ii ...)

  res <- resolve(res)
  res <- as.list(res)

  invisible(res)
}) # report()
