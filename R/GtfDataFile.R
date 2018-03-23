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
  extend(TabularTextFile(..., columnNames=columnNames), c("GtfDataFile", uses("SequenceContigsInterface")))
})

setMethodS3("as.character", "GtfDataFile", function(x, ...) {
  s <- NextMethod("as.character")
  s <- c(s, getSeqGenericSummary(x, ...))
  s
})

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
  organism <- Arguments$getOrganism(organism);

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
  pathname <- do.call(findAnnotationData, args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    args$pattern <- sprintf("%s[.]lnk$", args$pattern)
    pathname <- do.call(findAnnotationData, args=args);
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
  # Argument 'organism':
  organism <- Arguments$getOrganism(organism);

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


setMethodS3("getSeqLengths", "GtfDataFile", function(this, unique=FALSE, onlyIfCached=FALSE, force=FALSE, ...) {
  uniqify <- function(lens, unique=FALSE) {
    # Cache in memory
    this$.seqLengths <- lens

    if (!unique || length(lens) < 1L) return(lens)
    mstr(1)
    names <- names(lens)
    mstr(names)
    dups <- duplicated(names)
    mstr(dups)
    if (any(dups)) lens <- lens[!dups]
    mstr(lens)
    lens
  } # uniqify()

  pathname <- getPathname(this)

  # (a) Check for cached results in memory
  lens <- this$.seqLengths
  if (!force && !is.null(lens)) return(uniqify(lens, unique=unique))

  # (b) Check for cached results on file
  dirs <- c("aroma.seq", getOrganism(this))
  key <- list(method="getSeqLengths", class=class(this), pathname=pathname)
  if (!force) {
    lens <- loadCache(key=key, dirs=dirs)
    if (!is.null(lens)) return(uniqify(lens, unique=unique))
  }

  # (c) Scan file too expensive?
  if (!force && onlyIfCached) return(NULL)

  # (d) Scan file?
  con <- gzfile(pathname, open="r")
  on.exit(close(con))
  names <- NULL
  while (length(bfr <- readLines(con, n=10e3)) > 0L) {
    bfr <- grep("^#", bfr, value=TRUE, invert=TRUE)
    bfr <- gsub("\t.*", "", bfr)
    names <- c(names, bfr)
  }
  lens <- rep(NA_integer_, times=length(names))
  names(lens) <- names

  # Cache
  saveCache(lens, key=key, dirs=dirs)

  lens
})
