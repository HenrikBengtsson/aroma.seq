###########################################################################/**
# @RdocDefault doBowtie2
# @alias doBowtie2.FastqDataSet
#
# @title "Read alignment using the Bowtie v2 aligner"
#
# \description{
#  @get "title" based on [1].
# }
#
# \usage{
#   @usage doBowtie2
#   @usage doBowtie2,FastqDataSet
# }
#
# \arguments{
#  \item{dataSet, df}{A @see "FastqDataSet".}
#  \item{reference}{A @see "FastaReferenceFile" or a @see "Bowtie2IndexSet" specifying the genome reference to align the FASTQ reads to.}
#  \item{...}{Additional arguments passed to @see "Bowtie2Alignment".}
#  \item{verbose}{See @see "Verbose".}
# }
#
# \value{
#   Returns a @see "BamDataSet".
# }
#
# \references{
#  [1] Bowtie2, John Hopkins University, 2013.
#      \url{http://bowtie-bio.sourceforge.net/bowtie2/}
# }
#
# @author "HB"
#
# \seealso{
#  For more details, see @see "Bowtie2Alignment".
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("doBowtie2", "FastqDataSet", function(dataSet, reference, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'reference':
  if (inherits(reference, "FastaReferenceFile")) {
  } else if (inherits(reference, "BwaIndexSet")) {
  } else {
    throw("Argument 'reference' should either be of class 'FastaReferenceFile' or 'BwaIndexSet': ", class(reference)[1L]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Bowtie2");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Checking requirements
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Bowtie2/Check requirements");
  stopifnot(isCapableOf(aroma.seq, "bwa"));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Bowtie2 'aln' with options '-n 2' and '-q 40'.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "Bowtie2/Alignment");

  # Retrieve Bowtie2 index set?
  if (inherits(reference, "Bowtie2IndexSet")) {
    is <- reference;
  } else if (inherits(reference, "FastaReferenceFile")) {
    verbose && enter(verbose, "Bowtie2/Alignment/Retrieving index set");
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

  alg <- Bowtie2Alignment(dataSet, indexSet=is, ...);
  verbose && print(verbose, alg);

  bs <- process(alg, verbose=verbose);
  verbose && print(verbose, bs);

  verbose && exit(verbose);


  verbose && exit(verbose);

  bs;
}) # doBowtie2()


setMethodS3("doBowtie2", "default", function(...) {
  throw("The \"default\" method is still not implemented. Please see help('doBowtie2').");
})


############################################################################
# HISTORY:
# 2013-08-22
# o Created.
############################################################################
