###########################################################################/**
# @RdocClass FastaReferenceFile
#
# @title "The FastaReferenceFile class"
#
# \description{
#  @classhierarchy
#
#  A FastaReferenceFile object represents a FASTA reference file.
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
# \seealso{
#   ...
# }
#*/###########################################################################
setConstructorS3("FastaReferenceFile", function(...) {
  extend(GenericDataFile(...), "FastaReferenceFile",
    .seqLengths=NULL
  );
})

setMethodS3("as.character", "FastaReferenceFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  n <- nbrOfSeqs(this);
  s <- c(s, sprintf("Total sequence length: %d", getTotalSeqLengths(this)));
  s <- c(s, sprintf("Number of sequences: %s", n));
  s <- c(s, sprintf("Sequence names: [%d] %s", n, hpaste(getSeqNames(this))));

  class(s) <- class;
  s;
})


setMethodS3("getSeqLengths", "FastaReferenceFile", function(this, force=FALSE, ...) {
  seqLengths <- this$.seqLengths;
  if (force || is.null(seqLengths)) {
    seqLengths <- readSeqLengths(this, ...);
    this$.seqLengths <- seqLengths;
  }
  seqLengths;
})

setMethodS3("getTotalSeqLengths", "FastaReferenceFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...);
  res <- sum(as.numeric(seqLengths));
  if (res < .Machine$integer.max) {
    res <- as.integer(res);
  }
  res;
})

setMethodS3("getSeqNames", "FastaReferenceFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...);
  names(seqLengths);
})

setMethodS3("nbrOfSeqs", "FastaReferenceFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...);
  length(seqLengths);
})


# \seealso{
#   Internally, \code{fasta.info()} of \pkg{Biostrings} is used.
# }
setMethodS3("readSeqLengths", "FastaReferenceFile", function(this, ...) {
  pathname <- getPathname(this);
  seqLengths <- Biostrings::fasta.info(pathname);
  seqLengths;
}, private=TRUE)



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# BWA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
###########################################################################/** 
# @RdocMethod buildBwaIndexSet
#
# @title "Builds a BWA index files set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to @see "bwaIndex".}
#  \item{skip}{If @TRUE, the index files are not rebuilt if already available.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "R.filesets::GenericDataFileSet" consisting of the BWA index files.
# }
#
# @author
#*/########################################################################### 
setMethodS3("buildBwaIndexSet", "FastaReferenceFile", function(this, ..., skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  findIndexFiles <- function(prefix, exts=c("amb", "ann", "bwt", "pac", "sa"), ...) {
    pathnames <- sprintf("%s.%s", prefix, exts);
    names(pathnames) <- exts;
    pathnames[!file.exists(pathnames)] <- NA;
    pathnames;
  } # findIndexFiles()

  getIndexFileSet <- function(prefix, ...) {
    pathnames <- findIndexFiles(prefix);
    if (any(is.na(pathnames))) return(NULL);
    dfList <- lapply(pathnames, FUN=GenericDataFile);
  } # getIndexFileSet()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'skip':
  skip <- Arguments$getLogical(skip);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  } 
 

  verbose && enter(verbose, "Building BWA index files");

  pathnameFA <- getPathname(this);
  verbose && cat(verbose, "FASTA reference file to be indexed: ", pathnameFA);

  # The index prefix
  prefix <- bwaIndexPrefix(pathnameFA, ...);
  verbose && cat(verbose, "Prefix for index files: ", prefix);

  # Locate existing index files
  res <- tryCatch({
    BwaIndexSet$byPrefix(prefix);
  }, error=function(ex) BwaIndexSet());

  # Nothing todo?
  if (skip && isComplete(res)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && exit(verbose);
    return(res);
  }

  res <- bwaIndex(pathnameFA, indexPrefix=prefix, ..., verbose=less(verbose, 5));

  if (res != 0L) {
    throw("Failed to build BWA index. Return code: ", res);
  }

  res <- BwaIndexSet$byPrefix(prefix);

  # Sanity check
  stopifnot(!is.null(res));

  verbose && exit(verbose);

  res;
}) # buildBwaIndexSet()



############################################################################
# HISTORY:
# 2012-09-24
# o Now getTotalSeqLengths() returns a numeric if the result cannot be
#   held in an integer.
# o Added buildBwaIndexSet().
# 2012-06-28
# o Created.
############################################################################
