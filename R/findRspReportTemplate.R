setMethodS3("findRspReportTemplate", "Object", function(this, tags=NULL, type="(html|md|tex)", ..., paths="reports,rsp/", firstOnly=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  tags <- Arguments$getTags(tags);

  # Argument 'firstOnly':
  firstOnly <- Arguments$getLogical(firstOnly);

  # Argument 'type':
  type <- Arguments$getRegularExpression(type);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Finding RSP report templates");

  # Paths to search
  paths <- sapply(paths, FUN=function(path) {
    parts <- strsplit(path, split="::", fixed=TRUE)[[1]];
    if (length(parts) > 1L) {
      path <- system.file(parts[2L], package=parts[1L], mustWork=FALSE);
    }
    Arguments$getReadablePath(path, mustExist=FALSE);
  })
  paths <- paths[sapply(paths, FUN=isDirectory)];

  # Nothing to do?
  if (length(paths) == 0L) {
    return(NULL);
  }

  # Filename pattern to search for
  className <- class(this)[1L];
  fullname <- paste(c(className, tags), collapse=",");
  pattern <- sprintf("^%s.*[.]%s[.]rsp$", fullname, type);
  verbose && cat(verbose, "Filename pattern: ", pattern);

  pathnames <- c();
  for (pp in seq_along(paths)) {
    path <- paths[pp];
    verbose && enter(verbose, sprintf("Path #%d ('%s') of %d", pp, path, length(paths)));

    pathnamesT <- list.files(path=path, pattern=pattern, all.files=TRUE, full.names=TRUE);
    if (length(pathnamesT) == 0L) {
      verbose && exit(verbose);
      next;
    }

    if (firstOnly) {
      pathname <- pathnamesT[1L];
      verbose && exit(verbose);
      verbose && exit(verbose);
      return(pathname);
    }

    pathnames <- c(pathnames, pathnamesT);
    verbose && exit(verbose);
  } # for (path ...)


  verbose && cat(verbose, "Number of templates found: ", length(pathnames));
  verbose && print(verbose, pathnames);

  verbose && exit(verbose);

  pathnames;
}) # findRspReportTemplate()


setMethodS3("findRspReportTemplate", "FastqDataFile", function(this, ..., flavor="qrqc", paths=c("reports,rsp", "aroma.seq::reports,rsp")) {
  NextMethod("findRspReportTemplate", paths=paths);
}, protected=TRUE)

setMethodS3("findRspReportTemplate", "FastqDataSet", function(this, ..., flavor="qrqc", paths=c("reports,rsp", "aroma.seq::reports,rsp")) {
  NextMethod("findRspReportTemplate", paths=paths);
}, protected=TRUE)



############################################################################
# HISTORY:
# 2014-02-27
# o Added findRspReportTemplate() for FastqDataSet.
# 2013-11-12
# o Now findRspReportTemplate() expands 'paths' with format "pkg::path/to"
#   to system.file("path/to", package="pkg").
# o Added findRspReportTemplate() for FastqDataFile.
# 2012-12-06
# o Added findRspReportTemplate() for Object.
# o Created.
############################################################################
