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
setConstructorS3("GtfDataFile", function(...) {
  extend(TabularTextFile(...), "GtfDataFile")
})

setMethodS3("getOrganism", "GtfDataFile", function(this, ...) {
  path <- getPath(this);
  organism <- basename(path);
  organism;
})

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


  args <- list(pattern=pattern);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/organisms/<organism>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the fullname
  fullname <- paste(c(organism, tags), collapse=",");

  # Extract the name and the tags
  parts <- unlist(strsplit(fullname, split=",", fixed=TRUE));
  organism <- parts[1L];
  tags <- parts[-1L];

  # Search for "organisms/<organism>/.*[.]gtf$" files
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
  pathname <- findByOrganism(static, organism, ...);
  if (length(pathname) == 0)
    throw("Failed to located GTF file for organism: ", organism);
  res <- newInstance(static, pathname);
  organismR <- getOrganism(res);
  if (organismR != organism) {
    throw(sprintf("The located %s (%s) specifies an organism different from the requested one: %s != %s", class(res)[1L], getPathname(res), sQuote(organismR), sQuote(organism)));
  }
  res;
}, static=TRUE) # byOrganism()


############################################################################
# HISTORY:
# 2014-01-18
# o Created.
############################################################################
