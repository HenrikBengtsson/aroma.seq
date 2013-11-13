setMethodS3("report", "FastqDataFile", function(this, dataSet, ..., flavor="qrqc", type="md", outPath=".", verbose=FALSE) {
  require("R.rsp") || throw("Package not loaded: R.rsp");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'outPath':
  outPath <- Arguments$getWritablePath(outPath);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Generating report");
  verbose && cat(verbose, "Flavor: ", flavor);
  verbose && cat(verbose, "Data set: ", dataSet);
  verbose && print(verbose, this);

  # Locate report template
  pathname <- findRspReportTemplate(this, tags=flavor, type=type);
  if (length(pathname) == 0L) {
    throw("Failed to locate report template for this flavor and type: ", paste(c(sQuote(flavor), sQuote(type)), collapse=" & "));
  }
  verbose && cat(verbose, "Report template (source): ", pathname);

  # Output path
  outPath <- file.path(outPath, getFullName(this));
  outPath <- Arguments$getWritablePath(outPath);
  verbose && cat(verbose, "Report path: ", outPath);

  # Rename template
  filename <- basename(pathname);
  filenameT <- gsub("^[^,]*", getFullName(this), filename);
  pathnameT <- file.path(outPath, filenameT);
  copyFile(pathname, pathnameT, overwrite=TRUE);
  verbose && cat(verbose, "Report template (renamed): ", pathnameT);

  # Generate PDF report
  rspArgs <- list(this, dataSet=dataSet);
  pathnameR <- rfile(pathnameT, ..., workdir=outPath, verbose=less(verbose, 5));

  verbose && cat(verbose, "Generated report: ", pathnameR);

  verbose && exit(verbose);

  invisible(pathnameR);
}) # report()


setMethodS3("report", "FastqDataSet", function(this, dataSet=getFullName(this), ..., flavor="qrqc", outPath=file.path("reports", fullname(dataSet, tags=flavor)), verbose=FALSE) {
  require("R.rsp") || throw("Package not loaded: R.rsp");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'flavor':
  flavor <- match.arg(flavor);

  # Argument 'outPath':
  outPath <- Arguments$getWritablePath(outPath);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Generating reports");

  resList <- dsApply(this, function(df, dataSet, flavor, outPath, ..., verbose=FALSE) {
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'df':
    df <- Arguments$getInstanceOf(df, "FastqDataFile");

    # Argument 'dataSet':
    dataSet <- Arguments$getCharacter(dataSet);

    # Argument 'outPath':
    outPath <- Arguments$getWritablePath(outPath);

    # Argument 'flavor':
    flavor <- Arguments$getCharacter(flavor);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }

    res <- report(df, dataSet=dataSet, flavor=flavor, outPath=outPath, ..., verbose=less(verbose));

    res;
  }, dataSet=dataSet, flavor=flavor, outPath=outPath, verbose=verbose)

  verbose && exit(verbose);

  invisible(resList);
}) # report()


############################################################################
# HISTORY:
# 2013-11-12
# o SPEEDUP: Now report() for FastqDataSet supports parallel processing.
# 2013-07-18
# o BUG FIX: Incorrectly named argument ('dataSetSet').
# 2012-12-06
# o Added report() for FastqDataFile.
# o Created.
############################################################################
