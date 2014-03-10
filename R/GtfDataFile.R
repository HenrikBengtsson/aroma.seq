###########################################################################/**
# @RdocClass GtfDataFile
#
# @title "The GtfDataFile class"
#
# \description{
#  @classhierarchy
#
#  A GtfDataFile object represents a Gene Transfer Format (GTF) file.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::TabularTextFile".}
#   \item{columnNames}{Passed to @see "R.filesets::TabularTextFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Compression}{
#  The package supports compressed GTF files.
# }
#
# @author "HB"
#
# \seealso{
#   ...
# }
#*/###########################################################################
setConstructorS3("GtfDataFile", function(..., columnNames=FALSE) {
  extend(TabularTextFile(..., columnNames=columnNames), "GtfDataFile")
})

setMethodS3("as.character", "GtfDataFile", function(this, ...) {
  s <- NextMethod("as.character");
  seqNames <- getSeqNames(this, unique=TRUE, onlyIfCached=TRUE);
  s <- c(s, sprintf("Unique sequence names: %s [%d]", hpaste(seqNames), length(seqNames)));
  s;
}, protected=TRUE)


setMethodS3("getOrganism", "GtfDataFile", function(this, ...) {
  path <- getPath(this);
  organism <- basename(path);
  organism;
})


###########################################################################/**
# @RdocMethod byOrganism
# @aliasmethod findByOrganism
#
# @title "Locates a GTF file by organism"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{organism}{A @character string specifying for which organism a
#    file should be retrieved.}
#  \item{tags}{(not used) A @character @vector.}
#  \item{prefix}{(optional) A @character string specifying an optional
#    regular expression prefix to be prepended to \code{pattern} when
#    searching for the file.}
#  \item{pattern}{A @character string specifying a regular expression for
#    the file to be located.}
#  \item{...}{Additional arguments passed to the constructor of
#    @see "GtfDataFile" when instantiating the object.}
# }
#
# \value{
#   Returns a @see "GtfDataFile".
# }
#
# \seealso{
#   @seeclass
# }
#
# @author
#*/###########################################################################
setMethodS3("findByOrganism", "GtfDataFile", function(static, organism, tags=NULL, prefix=NULL, pattern="[.]gtf(|[.]gz)$", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'organism':
  organism <- Arguments$getCharacter(organism);

  # Argument 'prefix':
  if (!is.null(prefix)) {
    prefix <- Arguments$getRegularExpression(prefix);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/organisms/<organism>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the fullname
  fullname <- paste(c(organism, tags), collapse=",");

  # Extract the name and the tags
  parts <- unlist(strsplit(fullname, split=",", fixed=TRUE));
  organism <- parts[1L];
  tags <- parts[-1L];

  # Search for "organisms/<organism>/<prefix>.*[.]gtf$" files
  patternS <- pattern;
  if (!is.null(prefix)) patternS <- sprintf("%s.*%s", prefix, patternS);
  args <- list(
    set="organisms",
    name=organism,
    pattern=patternS,
    ...
  );
  pathname <- do.call("findAnnotationData", args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    args$pattern <- sprintf("%s[.]lnk$", args$pattern)
    pathname <- do.call("findAnnotationData", args=args);
    if (!is.null(pathname)) {
      # ..and expand it
      pathname <- Arguments$getReadablePathname(pathname, mustExist=FALSE);
      if (!isFile(pathname))
        pathname <- NULL;
    }
  }

  pathname;
}, static=TRUE, protected=TRUE) # findByOrganism()


setMethodS3("byOrganism", "GtfDataFile", function(static, organism, ...) {
  # Locate GTF file
  pathname <- findByOrganism(static, organism, ...);
  if (length(pathname) == 0L)
    throw("Failed to located GTF file for organism: ", organism);

  # Allocate object
  res <- newInstance(static, pathname, ..., .onUnknownArgs="ignore");

  # Validate
  organismR <- getOrganism(res);
  if (organismR != organism) {
    throw(sprintf("The located %s (%s) specifies an organism different from the requested one: %s != %s", class(res)[1L], getPathname(res), sQuote(organismR), sQuote(organism)));
  }
  res;
}, static=TRUE) # byOrganism()


setMethodS3("getSeqNames", "GtfDataFile", function(this, unique=FALSE, onlyIfCached=FALSE, force=FALSE, ...) {
  seqNames <- this$.seqNames;
  if (force || is.null(seqNames)) {
    if (!onlyIfCached) {
      pathname <- getPathname(this);
      con <- gzfile(pathname, open="r");
      on.exit(close(con));
      seqNames <- NULL;
      while (length(bfr <- readLines(con, n=10e3)) > 0L) {
        bfr <- gsub("\t.*", "", bfr);
        seqNames <- c(seqNames, bfr);
      }
      this$.seqNames <- seqNames;
    }
  }

  if (unique && length(seqNames) > 1L) {
    seqNames <- unique(seqNames);
    seqNames <- sort(seqNames);
    o <- order(nchar(seqNames), order(seqNames));
    seqNames <- seqNames[o];
  }

  seqNames;
})


############################################################################
# HISTORY:
# 2014-03-10
# o Added getSeqNames() for GtfDataFile, which now as.character() reports
#   on, iff already parsed.
# 2014-02-25
# o Now static byOrganism() no longer passes '...'.
# 2014-01-25
# o DOCUMENTATION: Added help for GtfDataFile$byOrganism().
# o Now static byOrganism() passes '...' also to the constructor.
# o Now GtfDataFile by default assumes that the GTF file has no
#   column headers.
# 2014-01-18
# o Created.
############################################################################
