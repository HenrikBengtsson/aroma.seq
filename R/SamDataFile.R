###########################################################################/**
# @RdocClass SamDataFile
#
# @title "The abstract SamDataFile class"
#
# \description{
#  @classhierarchy
#
#  A SamDataFile object represents a SAM file.
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#
# \references{
#    ...
# }
#
# \seealso{
#   An object of this class is typically part of an 
#   @see "SamDataSet".
# }
#*/###########################################################################
setConstructorS3("SamDataFile", function(...) {
  extend(GenericDataFile(...), "SamDataFile");
})


setMethodS3("as.character", "SamDataFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  class(s) <- class;
  s;
})


setMethodS3("convertToBamDataFile", "SamDataFile", function(this, path=getPath(this), ..., skip=!overwrite, overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePath(path);

  # Argument 'skip':
  skip <- Arguments$getLogical(skip);

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 
 

  verbose && enter(verbose, "Converting SAM file to a BAM file");
 
  pathname <- getPathname(this);
  verbose && cat(verbose, "SAM pathname: ", pathname);

  fullname <- getFullName(this);
  filenameBAM <- sprintf("%s.bam", fullname);
  pathnameBAM <- file.path(path, filenameBAM);
  verbose && cat(verbose, "BAM pathname: ", pathnameBAM);

  # Nothing to do?
  if (skip && isFile(pathnameBAM)) {
    verbose && cat(verbose, "Already converted. Skipping.");
    res <- BamDataFile(pathnameBAM);
    verbose && exit(verbose);
    return(res);
  }

  # Asserts
  stopifnot(getAbsolutePath(pathnameBAM) != getAbsolutePath(pathname));
  pathnameBAM <- Arguments$getWritablePathname(pathnameBAM, mustNotExist=!overwrite);

  # Converting
  verbose && enter(verbose, "Converting using Rsamtools");
  require("Rsamtools") || throw("Package not loaded: Rsamtools");
  pathnameBAMx <- gsub("[.]bam$", "", pathnameBAM);
  verbose && cat(verbose, "BAM destination: ", pathnameBAMx);
  # NB: Rsamtools::asBam() already writes atomically.
  pathnameD <- asBam(pathname, destination=pathnameBAMx, overwrite=overwrite);
  verbose && exit(verbose);

  res <- BamDataFile(pathnameBAM);

  verbose && exit(verbose);

  res;
}) # convertToBamDataFile()


############################################################################
# HISTORY:
# 2012-09-25
# o Created from BamDataFile.R.
############################################################################
