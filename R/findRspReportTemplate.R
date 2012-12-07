setMethodS3("findRspReportTemplate", "Object", function(this, tags=NULL, type=c("tex"), ..., paths="reports,rsp/", firstOnly=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  tags <- Arguments$getTags(tags);

  # Argument 'type':
  type <- match.arg(type);

  # Argument 'firstOnly':
  firstOnly <- Arguments$getLogical(firstOnly);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }  


  verbose && enter(verbose, "Finding RSP report templates");  

  # Paths to search
  paths <- sapply(paths, FUN=Arguments$getReadablePath, mustExist=FALSE);
  paths <- paths[sapply(paths, FUN=isDirectory)];

  # Nothing to do?
  if (length(paths) == 0L) {
    return(NULL);
  }

  # Filename pattern to search for
  className <- class(this)[1];
  fullname <- paste(c(className, tags), collapse=",");
  pattern <- sprintf("^%s.*.%s.rsp$", fullname, type);
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



############################################################################
# HISTORY:
# 2012-12-06
# o Added findRspReportTemplate().
# o Created.
############################################################################
