setMethodS3("findRspReportTemplate", "Object", function(this, tags=NULL, flavor=NULL, type="(html|md|tex)", ..., paths="reports,rsp/", firstOnly=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'tags':
  tags <- Arguments$getTags(tags)

  # Argument 'firstOnly':
  firstOnly <- Arguments$getLogical(firstOnly)

  # Argument 'flavor':

  # Argument 'type':
  type <- Arguments$getRegularExpression(type)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Finding RSP report templates")

  # Paths to search
  paths <- sapply(paths, FUN=function(path) {
    parts <- strsplit(path, split="::", fixed=TRUE)[[1]]
    if (length(parts) > 1L) {
      path <- system.file(parts[2L], package=parts[1L], mustWork=FALSE)
    }
    Arguments$getReadablePath(path, mustExist=FALSE)
  })
  paths <- paths[sapply(paths, FUN=isDirectory)]

  # Nothing to do?
  if (length(paths) == 0L) {
    return(NULL)
  }

  if (is.null(flavor)) {
    flavorPattern <- ""
  } else {
    flavorPattern <- paste(c("", flavor), collapse=",")
  }

  # Filename patterns to search for
  patterns <- sapply(class(this), FUN=function(className) {
    fullname <- paste(c(className, tags), collapse=",")
    sprintf("^%s.*%s[.]%s[.]rsp$", fullname, flavorPattern, type)
  })
  verbose && cat(verbose, "Filename patterns:")
  verbose && print(verbose, patterns)

  pathnames <- c()
  for (pp in seq_along(paths)) {
    path <- paths[pp]
    verbose && enter(verbose, sprintf("Path #%d ('%s') of %d", pp, path, length(paths)))

    pathnamesP <- c()
    for (qq in seq_along(patterns)) {
      pattern <- patterns[qq]
      verbose && enter(verbose, sprintf("Pattern #%d ('%s') of %d", qq, pattern, length(patterns)))
      pathnamesQ <- list.files(path=path, pattern=pattern, all.files=TRUE, full.names=TRUE)
      # Nothing found?
      if (length(pathnamesQ) == 0L) {
        verbose && exit(verbose)
        next
      }
      pathnamesP <- c(pathnamesP, pathnamesQ)
      verbose && exit(verbose)
      if (firstOnly) break
    } # for (qq ...)

    # Nothing found?
    if (length(pathnamesP) == 0L) {
      verbose && exit(verbose)
    }

    # Found something!

    # Nothing more to be done?
    if (firstOnly) {
      pathname <- pathnamesP[1L]
      verbose && exit(verbose)
      verbose && exit(verbose)
      return(pathname)
    }

    pathnames <- c(pathnames, pathnamesP)
    verbose && exit(verbose)
  } # for (pp ...)


  verbose && cat(verbose, "Number of templates found: ", length(pathnames))
  verbose && print(verbose, pathnames)

  verbose && exit(verbose)

  pathnames
}) # findRspReportTemplate()


setMethodS3("findRspReportTemplate", "FastqDataFile", function(this, ..., flavor="qrqc", paths=c("reports,rsp", "aroma.seq::reports,rsp")) {
  NextMethod("findRspReportTemplate", paths=paths)
}, protected=TRUE)

setMethodS3("findRspReportTemplate", "FastqDataSet", function(this, ...) {
  aFile <- getOneFile(this)
  findRspReportTemplate(aFile, ...)
}, protected=TRUE)
