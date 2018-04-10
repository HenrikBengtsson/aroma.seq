###########################################################################/**
# @RdocDefault doBWA
# @alias doBWA.AromaSeqDataFileSet
#
# @title "Read alignment using the Burrows-Wheeler Transform aligner (BWA)"
#
# \description{
#  @get "title" based on [1].
# }
#
# \usage{
#   @usage doBWA
#   @usage doBWA,AromaSeqDataFileSet
# }
#
# \arguments{
#  \item{dataSet, df}{A @see "FastqDataSet" or a @see "BamDataSet".}
#  \item{reference}{A @see "FastaReferenceFile" or a @see "BwaIndexSet" specifying the genome reference to align the FASTQ reads to.}
#  \item{...}{Additional arguments passed to @see "BwaAlignment".}
#  \item{verbose}{See @see "Verbose".}
# }
#
# \value{
#   Returns a @see "BamDataSet".
# }
#
# \references{
#   [1] Li H. and Durbin R., \emph{Fast and accurate short read alignment
#       with Burrows-Wheeler Transform}. Bioinformatics, 2009.\cr
#   [2] Li H. and Durbin R., \emph{Fast and accurate long-read alignment
#       with Burrows-Wheeler Transform}. Bioinformatics, 2010.\cr
# }
#
# @author "HB"
#
# \seealso{
#  For more details, see @see "BwaAlignment".
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("doBWA", "AromaSeqDataFileSet", function(dataSet, reference, ..., verbose=FALSE) {
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


  verbose && enter(verbose, "BWA");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Checking requirements
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "BWA/Check requirements");
  stop_if_not(isCapableOf(aroma.seq, "bwa"));
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # BWA 'aln' with options '-n 2' and '-q 40'.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "BWA/Alignment");

  # Retrieve BWA index set?
  if (inherits(reference, "BwaIndexSet")) {
    is <- reference;
  } else if (inherits(reference, "FastaReferenceFile")) {
    verbose && enter(verbose, "BWA/Alignment/Retrieving index set");
    fa <- reference;
    verbose && print(verbose, fa);
    is <- buildBwaIndexSet(fa, verbose=verbose);
    verbose && print(verbose, is);
    verbose && exit(verbose);
    # Not needed anymore
    fa <- NULL;
  }
    # Not needed anymore
  reference <- NULL;

  alg <- BwaAlignment(dataSet, indexSet=is, ...);
  verbose && print(verbose, alg);

  bs <- process(alg, verbose=verbose);
  verbose && print(verbose, bs);

  verbose && exit(verbose);


  verbose && exit(verbose);

  bs;
}) # doBWA()


setMethodS3("doBWA", "default", function(...) {
  throw("The \"default\" method is still not implemented. Please see help('doBWA').")
})
