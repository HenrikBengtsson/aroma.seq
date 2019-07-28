###########################################################################/**
# @RdocClass SraDataSet
#
# @title "The SraDataSet class"
#
# \description{
#  @classhierarchy
#
#  An SraDataSet object represents a set of @see "SraDataFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "SraDataFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SraDataSet", function(files=NULL, ...) {
  extend(GenericDataFileSet(files=files, ...), c("SraDataSet", uses("AromaSeqDataFileSet")))
})


setMethodS3("getOrganism", "SraDataSet", function(this, depth=getDepth(this)-1L, ...) {
  path <- getPath(this)
  path <- getParent(path, depth=depth)
  organism <- basename(path)
  organism <- Arguments$getCharacter(organism, length=c(1L, 1L))
  organism
}, protected=TRUE)


setMethodS3("getDepth", "SraDataSet", function(this, ...) {
  path <- getPath(this, absolute=FALSE)
  parts <- unlist(strsplit(path, split="/", fixed=TRUE))
  nparts <- length(parts)
  depth <- nparts - 2L
  depth <- Arguments$getInteger(depth, range=c(0,Inf))
  depth
}, protected=TRUE)


setMethodS3("byPath", "SraDataSet", function(static, ..., pattern="[.](sra|SRA)$") {
  NextMethod("byPath", pattern=pattern)
}, static=TRUE)



setMethodS3("findByName", "SraDataSet", function(static, name, tags=NULL, organism=NULL, ..., paths="sraData", pattern="[.](sra|SRA)$") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'organism':
  if (!is.null(organism)) {
    organism <- Arguments$getOrganism(organism)
  }

  NextMethod("findByPath", subdirs=organism, paths=paths, pattern=pattern)
}, static=TRUE)



setMethodS3("byName", "SraDataSet", function(static, name, tags=NULL, organism=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'organism':
  if (!is.null(organism)) {
    organism <- Arguments$getOrganism(organism)
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Setting up ", class(static)[1L], " by name")

  verbose && cat(verbose, "Organism: ", organism)

  suppressWarnings({
    paths <- findByName(static, name, tags=tags, organism=organism,
                                                   firstOnly=FALSE, ...)
  })
  if (is.null(paths)) {
    path <- file.path(paste(c(name, tags), collapse=","), organism)
    throw("Cannot create ", class(static)[1], ".  No such directory: ", path)
  }

  verbose && cat(verbose, "Paths to possible data sets:")
  verbose && print(verbose, paths)

  # Record all exception
  exList <- list()

  res <- NULL
  for (kk in seq_along(paths)) {
    path <- paths[kk]
    verbose && enter(verbose, sprintf("Trying path #%d of %d", kk, length(paths)))
    verbose && cat(verbose, "Path: ", path)

    tryCatch({
      suppressWarnings({
        res <- byPath(static, path=path, ..., verbose=verbose)
      })
    }, error = function(ex) {
      exList <<- append(exList, list(list(exception=ex, path=path)))

      verbose && cat(verbose, "Data set could not be setup for this path, because:")
      verbose && cat(verbose, ex$message)
    })

    if (!is.null(res)) {
      if (length(res) > 0) {
        verbose && cat(verbose, "Successful setup of data set.")
        verbose && exit(verbose)
        break
      }
    }

    verbose && exit(verbose)
  } # for (kk ...)

  if (is.null(res)) {
    exMsgs <- sapply(exList, FUN=function(ex) {
      sprintf("%s (while trying '%s').",
                   ex$exception$message, ex$path)
    })
    exMsgs <- sprintf("(%d) %s", seq_along(exMsgs), exMsgs)
    exMsgs <- paste(exMsgs, collapse="  ")
    msg <- sprintf("Failed to setup a data set for any of %d data directories located. The following reasons were reported: %s", length(paths), exMsgs)
    verbose && cat(verbose, msg)
    throw(msg)
  }

  verbose && exit(verbose)

  res
}, static=TRUE)



setMethodS3("fastqDump", "SraDataSet", function(this, path=".", ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePath(path)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Writing SRA set as FASTQ set")
  verbose && cat(verbose, "Input set:")
  verbose && print(verbose, this)
  verbose && cat(verbose, "Output path: ", path)

  fqs <- lapply(this, FUN=fastqDump, path=path, ..., verbose=verbose)
  fqs <- FastqDataSet(fqs)
  verbose && print(verbose, fqs)

  verbose && exit(verbose)

  fqs
})
