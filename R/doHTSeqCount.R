###########################################################################/**
# @RdocDefault doHTSeqCount
# @alias doHTSeqCount.BamDataSet
#
# @title "Counting reads in features"
#
# \description{
#  @get "title" based on [1].
# }
#
# \usage{
#   @usage doHTSeqCount
#   @usage doHTSeqCount,BamDataSet
# }
#
# \arguments{
#  \item{dataSet}{A @see "BamDataSet".}
#  \item{transcripts}{A @see "GtfDataFile".}
#  \item{...}{Additional arguments passed to @see "HTSeqCounting".}
#  \item{verbose}{See @see "Verbose".}
# }
#
# \value{
#   Returns a @see "HTSeqCountDataSet".
# }
#
# \references{
#  [1] Simon Anders, \emph{HTSeq: Analysing high-throughput sequencing
#      data with Python}, EMBL, Jan 2014.
#      \url{http://www-huber.embl.de/users/anders/HTSeq/} \cr 
# }
#
# @author "HB"
#
# \seealso{
#  For more details, see @see "HTSeqCounting".
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("doHTSeqCount", "BamDataSet", function(dataSet, transcripts, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getInstanceOf(dataSet, "BamDataSet");

  # Argument 'transcripts':
  transcripts <- Arguments$getInstanceOf(transcripts, "GtfDataFile"); 

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "HTSeqCounting");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Checking requirements
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "HTSeqCounting/Check requirements");
  stop_if_not(isCapableOf(aroma.seq, "htseq"));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # HTSeqCounting
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "HTSeqCounting/Alignment");
  htc <- HTSeqCounting(dataSet, transcripts=transcripts, ...);
  verbose && print(verbose, htc);

  counts <- process(htc, verbose=verbose);
  verbose && print(verbose, counts);

  verbose && exit(verbose);


  verbose && exit(verbose);

  counts;
}) # doHTSeqCount()


setMethodS3("doHTSeqCount", "default", function(...) {
  throw("In order to use doHTSeqCount(), the data need to be aligned into a BamDataFileSet first.");
})
