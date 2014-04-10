###########################################################################/**
# @RdocDefault doTopHat2
# @alias doTopHat2.FastqDataSet
#
# @title "Read alignment using the TopHat v2 aligner"
#
# \description{
#  @get "title" based on [1].
# }
#
# \usage{
#   @usage doTopHat2
#   @usage doTopHat2,FastqDataSet
# }
#
# \arguments{
#  \item{dataSet, df}{A @see "FastqDataSet".}
#  \item{reference}{A @see "FastaReferenceFile" or a @see "Bowtie2IndexSet" specifying the genome reference to align the FASTQ reads to.}
#  \item{transcripts}{A @see "GtfDataFile" or @NULL specifying known transcripts and/or gene model annotations.}
#  \item{...}{Additional arguments passed to @see "TopHat2Alignment".}
#  \item{verbose}{See @see "Verbose".}
# }
#
# \value{
#   Returns a @see "BamDataSet".
# }
#
# \references{
#  [1] TopHat, University of Maryland, 2013.
#      \url{http://http://tophat.cbcb.umd.edu/} \cr
#  [2] Trapnell et al. \emph{Differential gene and transcript expression
#      analysis of RNA-seq experiments with TopHat and Cufflinks}.
#      Nat Protoc, 2012.\cr
# }
#
# @author "HB"
#
# \seealso{
#  For more details, see @see "TopHat2Alignment".
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("doTopHat2", "FastqDataSet", function(dataSet, reference, transcripts, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'reference':
  if (inherits(reference, "FastaReferenceFile")) {
  } else if (inherits(reference, "Bowtie2IndexSet")) {
  } else {
    throw("Argument 'reference' should either be of class 'FastaReferenceFile' or 'BwaIndexSet': ", class(reference)[1L]);
  }

  # Argument 'reference':
  if (is.null(transcripts)) {
  } else if (inherits(transcripts, "GtfDataFile")) {
  } else {
    throw("Argument 'transcripts' should be either of class 'GtfDataFile' or explicitly NULL': ", class(transcripts)[1L]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "TopHat2");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Checking requirements
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "TopHat2/Check requirements");
  stopifnot(isCapableOf(aroma.seq, "tophat2"));
  stopifnot(isCapableOf(aroma.seq, "bowtie2"));
  stopifnot(isCapableOf(aroma.seq, "samtools"));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # TopHat2
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "TopHat2/Alignment");

  # Retrieve Bowtie2 index set?
  if (inherits(reference, "Bowtie2IndexSet")) {
    is <- reference;
  } else if (inherits(reference, "FastaReferenceFile")) {
    verbose && enter(verbose, "TopHat2/Alignment/Retrieving index set");
    fa <- reference;
    verbose && print(verbose, fa);
    is <- buildBowtie2IndexSet(fa, verbose=verbose);
    verbose && print(verbose, is);
    verbose && exit(verbose);
    # Not needed anymore
    fa <- NULL;
  }
    # Not needed anymore
  reference <- NULL;

  alg <- TopHat2Alignment(dataSet, indexSet=is, transcripts=transcripts, ...);
  verbose && print(verbose, alg);

  bams <- process(alg, verbose=verbose);
  verbose && print(verbose, bams);

  verbose && exit(verbose);


  verbose && exit(verbose);

  bams;
}) # doTopHat2()


setMethodS3("doTopHat2", "default", function(...) {
  throw("The \"default\" method is still not implemented. Please see help('doTopHat2').");
})


############################################################################
# HISTORY:
# 2014-04-10
# o Now doTopHat2() requires that 'transcripts' is explicitly specified.
#   This is to protect against the mistake when the user forgets to
#   specify the transcripts.
# 2013-11-02
# o Created from doBowtie2().
############################################################################
