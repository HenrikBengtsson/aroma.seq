###########################################################################/**
# @RdocGeneric readSeqFile
# @alias readSeqFile.FastqDataFile
# @alias readFASTQSummary
# @alias readFASTQSummary.FastqDataFile
#
# @title "Reads and summarizes a FASTQ file"
#
# \description{
#   @get "title"
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "qrqc::readSeqFile".}
#  \item{seed}{An (optional) @integer specifying the random seed to be
#     set before sampling indices.  The random seed is set to its original
#     state when exiting.  If @NULL, it is not set.}
#  \item{cache}{If @TRUE, memoization is used, otherwise not.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "qrqc::FASTQSummary-class" object.
# }
#
# \seealso{
#   Internally @see "qrqc::readSeqFile" is utilized.
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("readSeqFile", "FastqDataFile", function(this, ..., seed=NULL, cache=TRUE, verbose=FALSE) {
  # NOTE: Do not attach 'qrqc', because then qrqc::readSeqFile()
  # will override generic aroma.seq::readSeqFile().
  readSeqFile <- qrqc::readSeqFile;


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'seed':
  if (!is.null(seed)) {
    seed <- Arguments$getInteger(seed);
  }

  # Argument 'cache':
  cache <- Arguments$getLogical(cache);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Reading FASTQSummary");

  pathname <- getPathname(this);
  verbose && print(verbose, this);
  pathname <- Arguments$getReadablePathname(pathname);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for cached results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    key <- list(method="readSeqFile", class=class(this)[1L],
                pathname=pathname, ..., seed=seed);
    dirs <- c("aroma.seq", "readSeqFile");
    data <- loadCache(key=key, dirs=dirs);
    if (!is.null(data)) {
      verbose && cat(verbose, "Found cached results. Skipping.");
      verbose && exit(verbose);
      return(data);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Set the random seed
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (!is.null(seed)) {
    verbose && enter(verbose, "Setting (temporary) random seed");
    oldRandomSeed <- NULL;
    if (exists(".Random.seed", mode="integer")) {
      oldRandomSeed <- get(".Random.seed", mode="integer");
    }
    on.exit({
      if (!is.null(oldRandomSeed)) {
        .Random.seed <<- oldRandomSeed;
      }
    }, add=TRUE);
    verbose && cat(verbose, "The random seed will be reset to its original state afterward.");
    verbose && cat(verbose, "Seed: ", seed);
    set.seed(seed);
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "qrqc::readSeqFile()");
  data <- readSeqFile(pathname, type="fastq", ...);
  verbose && print(verbose, data);
  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Cache results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    saveCache(data, key=key, dirs=dirs);
    verbose && cat(verbose, "Saved results to cache.");
  }

  verbose && exit(verbose);

  data;
}) # readSeqFile()


setMethodS3("readFASTQSummary", "FastqDataFile", function(this, ...) {
  readFASTQSummary(...);
}, protected=TRUE)



## setOldClass("FastqDataFile")
