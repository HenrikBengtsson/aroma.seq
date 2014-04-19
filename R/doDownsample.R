###########################################################################/**
# @RdocDefault doDownsample
# @alias doDownsample.BamDataSet
# @alias doDownsample.FastqDataSet
#
# @title "Generates a downsampled FASTQ or BAM data set"
#
# \description{
#  @get "title".
# }
#
# \usage{
#   @usage doDownsample
#   @usage doDownsample,BamDataSet
#   @usage doDownsample,FastqDataSet
# }
#
# \arguments{
#  \item{dataSet}{A @see "BamDataSet" or @see "FastqDataSet".}
#  \item{subset}{An @integer specifying the total number of reads to sample,
#    or a @double specifying the fraction of total number of reads to sample.}
#  \item{...}{Additional arguments passed to specific downsampler, e.g.
#     @see "BamDownsampler" and @see "FastqDownsampler".}
#  \item{verbose}{See @see "Verbose".}
# }
#
# \value{
#   Returns a @see "GenericDataFileSet" of the same class
#   as the input data set \code{dataSet}.
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("doDownsample", "BamDataSet", function(dataSet, subset=1e6, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Downsampling");
  ds <- BamDownsampler(dataSet, subset=subset, ...);
  verbose && print(verbose, ds);
  dsOut <- process(ds, verbose=verbose);
  verbose && print(verbose, dsOut);
  verbose && exit(verbose);

  dsOut;
}) # doDownsample()


setMethodS3("doDownsample", "FastqDataSet", function(dataSet, subset=subset, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Downsampling");
  ds <- FastqDownsampler(dataSet, ...);
  verbose && print(verbose, ds);
  dsOut <- process(ds, verbose=verbose);
  verbose && print(verbose, dsOut);
  verbose && exit(verbose);

  dsOut;
}) # doDownsample()


setMethodS3("doDownsample", "default", function(...) {
  throw("The \"default\" method is still not implemented. Please see help('doDownsample').");
})


############################################################################
# HISTORY:
# 2014-04-18
# o Added doDownsample() for BamDataSet and FastqDataSet.
# o Created.
############################################################################
