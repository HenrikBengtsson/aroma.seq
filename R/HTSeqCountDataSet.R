###########################################################################/**
# @RdocClass HTSeqCountDataSet
#
# @title "The HTSeqCountDataSet class"
#
# \description{
#  @classhierarchy
#
#  An HTSeqCountDataSet object represents a set of @see "HTSeqCountDataFile":s.
# }
#
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "HTSeqCountDataFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::TabularTextFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \references{
#  [1] Simon Anders, \emph{HTSeq: Analysing high-throughput sequencing
#      data with Python}, EMBL, Jan 2014.
#      \url{http://www-huber.embl.de/users/anders/HTSeq/} \cr
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("HTSeqCountDataSet", function(files=NULL, ...) {
  extend(TabularTextFileSet(files=files, ...), c("HTSeqCountDataSet", uses("AromaSeqDataFileSet")));
})


setMethodS3("getDepth", "HTSeqCountDataSet", function(this, ...) {
  1L;
}, protected=TRUE)


setMethodS3("getOrganism", "HTSeqCountDataSet", function(this, ...) {
  organism <- directoryItem(this, "organism");
  organism <- Arguments$getCharacter(organism, length=c(1L, 1L));
  organism;
}, protected=TRUE)


setMethodS3("byPath", "HTSeqCountDataSet", function(static, ..., recursive=FALSE, pattern="[.](count)$", verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'recursive':
  recursive <- Arguments$getLogical(recursive);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Setting up ", class(static)[1L], " by path");
  verbose && cat(verbose, "Recursive: ", recursive);
  verbose && cat(verbose, "Filename pattern: ", pattern);

  res <- NextMethod("byPath", recursive=recursive, pattern=pattern);

  verbose && exit(verbose);

  res;
}, protected=TRUE)


setMethodS3("findByName", "HTSeqCountDataSet", function(static, name, tags=NULL, organism=NULL, ..., paths="htseqCountData", pattern="[.](count)$") {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'organism':
  if (!is.null(organism)) {
    organism <- Arguments$getCharacter(organism);
  }

  NextMethod("findByPath", subdirs=organism, paths=paths, pattern=pattern);
}, static=TRUE)



setMethodS3("byName", "HTSeqCountDataSet", function(static, name, tags=NULL, organism=NULL, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'organism':
  if (!is.null(organism)) {
    organism <- Arguments$getCharacter(organism);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Setting up ", class(static)[1L], " by name");

  verbose && cat(verbose, "Organism: ", organism);

  suppressWarnings({
    paths <- findByName(static, name, tags=tags, organism=organism,
                                      firstOnly=FALSE, ...);
  })
  if (is.null(paths)) {
    path <- file.path(paste(c(name, tags), collapse=","), organism);
    throw("Cannot create ", class(static)[1], ".  No such directory: ", path);
  }

  verbose && cat(verbose, "Paths to possible data sets:");
  verbose && print(verbose, paths);

  # Record all exception
  exList <- list();

  res <- NULL;
  for (kk in seq_along(paths)) {
    path <- paths[kk];
    verbose && enter(verbose, sprintf("Trying path #%d of %d", kk, length(paths)));
    verbose && cat(verbose, "Path: ", path);

    tryCatch({
      suppressWarnings({
        res <- byPath(static, path=path, ..., verbose=verbose);
      });
    }, error = function(ex) {
      exList <<- append(exList, list(list(exception=ex, path=path)));

      verbose && cat(verbose, "Data set could not be setup for this path, because:");
      verbose && cat(verbose, ex$message);
    });

    if (!is.null(res)) {
      if (length(res) > 0) {
        verbose && cat(verbose, "Successful setup of data set.");
        verbose && exit(verbose);
        break;
      }
    }

    verbose && exit(verbose);
  } # for (kk ...)

  if (is.null(res)) {
    exMsgs <- sapply(exList, FUN=function(ex) {
      sprintf("%s (while trying '%s').",
                   ex$exception$message, ex$path);
    });
    exMsgs <- sprintf("(%d) %s", seq_along(exMsgs), exMsgs);
    exMsgs <- paste(exMsgs, collapse="  ");
    msg <- sprintf("Failed to setup a data set for any of %d data directories located. The following reasons were reported: %s", length(paths), exMsgs);
    verbose && cat(verbose, msg);
    throw(msg);
  }

  verbose && exit(verbose);

  res;
}, static=TRUE)



###########################################################################/**
# @RdocMethod readDGE
# @alias readDGE
# @aliasmethod extractDGEList
# @alias extractDGEList
#
# @title "Reads all digital gene expression (DGE) data"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{labels}{A @character string specifying sample names.}
#  \item{...}{Additional arguments passed to @see "edgeR::readDGE".}
# }
#
# \value{
#   Returns a @see "edgeR::DGEList" object.
# }
#
# @author "HB"
#*/###########################################################################
setMethodS3("readDGE", "HTSeqCountDataSet", function(this, labels=getFullNames(this), ...) {
  use("edgeR")
  ns <- getNamespace("edgeR");
  readDGE <- get("readDGE", envir=ns, mode="function");

  pathnames <- getPathnames(this);

  # Argument 'labels':
  labels <- Arguments$getCharacters(labels, length=rep(length(pathnames), times=2L));

  readDGE(pathnames, labels=labels, ...);
}, protected=TRUE)


setMethodS3("extractDGEList", "HTSeqCountDataSet", function(this, ...) {
  readDGE(this, ...);
})



############################################################################
# HISTORY:
# 2014-01-24
# o Created.
############################################################################
