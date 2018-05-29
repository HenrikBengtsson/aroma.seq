###########################################################################/**
# @RdocDefault systemJava
#
# @title "Calls the Java executable"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Additional arguments specifying Java command line switches.}
#   \item{.fmtArg}{A named @character @vector specifying how command line
#    options are formatted, e.g. \code{"-f bar"} vs  \code{"--foo=bar"}.}
#   \item{.fake}{If @TRUE, the executable is not called.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemJava", "default", function(..., jvmArgs=Sys.getenv("JVM_ARGS"), .fmtArg=c("( |=)$"="%s%s", ".*"="%s %s"), .fake=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments '...':
  args <- list(...)

  # Arguments 'jvmArgs':
  if (is.character(jvmArgs)) {
    jvmArgs <- trim(jvmArgs)
    jvmArgs <- jvmArgs[nzchar(jvmArgs)]
  }

  # Argument 'autoDash':
  .fmtArg <- Arguments$getCharacters(.fmtArg, useNames=TRUE)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Calling Java executable")

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  bin <- findJava(verbose=less(verbose, 50))
  verbose && cat(verbose, "Java executable:")
  verbose && print(verbose, bin)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Split up '...' arguments by system2() and Java executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  keys <- names(args)
  if (is.null(keys)) {
    system2Args <- NULL
  } else {
    keep <- is.element(keys, names(formals(base::system2)))
    system2Args <- args[keep]
    args <- args[!keep]
  }

  verbose && cat(verbose, "Arguments passed to system2():")
  verbose && str(verbose, system2Args)

  verbose && cat(verbose, "Arguments passed to java:")
  verbose && str(verbose, args)

  verbose && cat(verbose, "JVM arguments passed to java:")
  verbose && str(verbose, jvmArgs)

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Adjust command line switches
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Adjusting command-line options")

  if (verbose) verbose <- less(verbose, 10)
  verbose && enter(verbose, "Adjusting patterns")
  # Default patterns
  defFmtArg <- c("( |=)$"="%s%s", ".*"="%s %s")

  # Insert default ones at the beginning, unless specified.
  if (!is.element("*", .fmtArg)) {
    .fmtArg <- c("*", .fmtArg)
  }

  at <- match("*", .fmtArg)[1L]
  .fmtArg <- R.utils::insert(.fmtArg, ats=at+1L, defFmtArg)[-at]
  .fmtArg <- rev(.fmtArg)
  dups <- duplicated(names(.fmtArg))
  .fmtArg <- .fmtArg[!dups]
  .fmtArg <- rev(.fmtArg)

  verbose && cat(verbose, "Patterns:")
  verbose && print(verbose, .fmtArg)
  verbose && exit(verbose)
  if (verbose) verbose <- more(verbose, 10)

  if (verbose) verbose <- less(verbose, 10)
  patterns <- names(.fmtArg)
  options <- NULL
  for (kk in seq_along(patterns)) {
    patternKK <- patterns[kk]
    fmtArgKK <- .fmtArg[kk]
    keys <- names(args)

    idxs <- grep(patternKK, keys)
    # No match?
    if (length(idxs) == 0L)
      next

    keysKK <- keys[idxs]
    argsKK <- args[idxs]

    verbose && cat(verbose, "Pattern: ", patternKK)
    verbose && cat(verbose, "Format: ", fmtArgKK)
    verbose && cat(verbose, "Keys: ", hpaste(keysKK))
    verbose && cat(verbose, "Values: ", hpaste(argsKK))
    optionsKK <- sprintf(fmtArgKK, keysKK, argsKK)
    # Undo non-named arguments
    hasNoName <- (nchar(keysKK) == 0L)
    optionsKK[hasNoName] <- argsKK[hasNoName]
    verbose && cat(verbose, "Options: ", hpaste(optionsKK))

    options <- c(options, optionsKK)

    # Done?
    args <- args[-idxs]
    if (length(args) == 0L)
      break
  } # for (kk ...)
  if (verbose) verbose <- more(verbose, 10)

  # Prepend JVM_ARGS
  options <- c(jvmArgs, options)

  options <- trim(options)
  options <- options[nchar(options) > 0L]
  verbose && cat(verbose, "Command line options:")
  verbose && print(verbose, options)

  verbose && exit(verbose)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # System call
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && cat(verbose, "Java system call:")
  cmd <- sprintf("%s %s", bin, paste(options, collapse=" "))
  verbose && print(verbose, cmd)
  verbose && str(verbose, system2Args)

  verbose && enter(verbose, "system2() call")
  callArgs <- list(command=bin, args=options)
  callArgs <- c(callArgs, system2Args)
  verbose && str(verbose, callArgs)
  if (!.fake) {
    res <- do.call(base::system2, callArgs)
  } else {
    res <- "<fake run>"
  }
  verbose && exit(verbose)

  verbose && exit(verbose)

  res
}) # systemJava()


###########################################################################/**
# @RdocDefault systemJavaJar
#
# @title "Calls the Java jar executable"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathname}{The JAR file.}
#   \item{...}{Additional arguments specifying Java command line switches.}
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("systemJavaJar", "default", function(pathname, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'pathname':
  pathname <- Arguments$getReadablePathname(pathname)

  # Quote jar pathname just in case
  systemJava("-jar "=sprintf("\"%s\"", pathname), ...)
})
