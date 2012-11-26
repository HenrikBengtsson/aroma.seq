###########################################################################/**
# @RdocClass SamDataSet
#
# @title "The SamDataSet class"
#
# \description{
#  @classhierarchy
#
#  An SamDataSet object represents a set of @see "SamDataFile":s.
# }
# 
# @synopsis
#
# \arguments{
#   \item{files}{A @list of @see "SamDataFile":s.}
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFileSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
# 
# @author
#*/###########################################################################
setConstructorS3("SamDataSet", function(files=NULL, ...) {
  extend(GenericDataFileSet(files=files, ...), "SamDataSet");
})


setMethodS3("validate", "SamDataSet", function(this, ...) {
  NextMethod("validate");
}, protected=TRUE)


setMethodS3("getDepth", "SamDataSet", function(this, ...) {
  1L;
}, protected=TRUE);


setMethodS3("byPath", "SamDataSet", function(static, ..., pattern="[.](sam|SAM)$") {
  NextMethod("byPath", pattern=pattern);
}, static=TRUE)



###########################################################################/** 
# @RdocMethod convertToBamDataSet
#
# @title "Converts the SAM set into a BAM set"
#
# \description{
#   @get "title", where each BAM file is sorted and indexed.
# }
#
# @synopsis
#
# \arguments{
#  \item{path}{The destination path.}
#  \item{...}{Additional arguments passed to \code{convertToBamDataFile()} 
#   for each @see "SamDataFile".}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "BamDataSet".
# }
#
# \seealso{
#   Internally @see "Rsamtools::asBam" is utilized.
# }
#
# @author
#*/###########################################################################  
setMethodS3("convertToBamDataSet", "SamDataSet", function(this, path=getPath(this), ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePath(path);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 
 

  verbose && enter(verbose, "Converting SAM set to a BAM set");
 
  verbose && cat(verbose, "SAM data set:");
  verbose && print(verbose, this);
  verbose && cat(verbose, "BAM path: ", path);

  for (ii in seq_along(this)) {
    sf <- getFile(this, ii);
    verbose && enter(verbose, sprintf("File #%d ('%s') of %d", ii, getName(sf), length(this)));
    bf <- convertToBamDataFile(sf, path=path, ..., verbose=less(verbose,1));
    verbose && exit(verbose);
  } # for (ii ...)

  bs <- BamDataSet$byPath(path);
  bs <- extract(bs, getFullNames(this), onMissing="error");
  verbose && print(verbose, bs);

  ## TODO: Assert completeness

  verbose && exit(verbose);

  bs;
}) # convertToBamDataSet()


############################################################################
# HISTORY:
# 2012-11-26
# o ROBUSTNESS: Now convertToBamDataSet() for SamDataSet asserts that
#   the output set is complete.
# 2012-09-25
# o Created from BamDataSet.R.
############################################################################
