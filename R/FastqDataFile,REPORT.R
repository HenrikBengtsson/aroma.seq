setMethodS3("report", "FastqDataFile", function(this, dataSetSet, ..., flavor="qrqc", outPath=".", verbose=FALSE) {
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

  paths <- c("reports,rsp", system.file("reports,rsp", package="aroma.seq"));
  pathname <- findRspReportTemplate(this, tags=flavor, paths=paths);
  if (length(pathname) == 0L) {
    throw("Failed to locate report template.");
  }
  verbose && cat(verbose, "Report template (source): ", pathname);

  # Rename template
  filename <- basename(pathname);
  filenameT <- gsub("^[^,]*", getFullName(this), filename);
  pathnameT <- file.path(outPath, filenameT);
  copyFile(pathname, pathnameT, overwrite=TRUE);
  verbose && cat(verbose, "Report template (renamed): ", pathnameT);

  # WORKAROUND: Force absolute path
  thisT <- newInstance(this, getAbsolutePath(getPathname(this)));
  rspArgs <- list(thisT, dataSet=dataSet);

  # Generate PDF report
  pathnameR <- rsp(pathnameT, ..., outPath=outPath, verbose=less(verbose, 5));

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

  resList <- vector("list", length=length(this));
  for (ii in seq_along(this)) {
    df <- getFile(this, ii);
    name <- getFullName(df);
    verbose && enter(verbose, sprintf("File #%d ('%s') of %d", ii, name, length(this)));
   
    resList[[ii]] <- report(df, dataSet=dataSet, flavor=flavor, outPath=outPath, ..., verbose=less(verbose));
    verbose && exit(verbose);
  } # for (ii ...)

  verbose && exit(verbose);

  invisible(resList);
}) # report()


############################################################################
# HISTORY:
# 2012-12-06
# o Added report() for FastqDataFile.
# o Created.
############################################################################
