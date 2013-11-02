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
# @author "HB"
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

  s <- NextMethod("as.character");
  class <- class(s);

  n <- nbrOfSeqs(this, onlyIfCached=TRUE);
  s <- c(s, sprintf("Total sequence length: %.0f", getTotalSeqLengths(this, onlyIfCached=TRUE)));
  s <- c(s, sprintf("Number of sequences: %d", n));
  s <- c(s, sprintf("Sequence names: [%d] %s", n, hpaste(getSeqNames(this, onlyIfCached=TRUE))));

  class(s) <- class;
  s;
}, protected=TRUE)


setMethodS3("getSeqLengths", "FastaReferenceFile", function(this, force=FALSE, onlyIfCached=FALSE, ...) {
  seqLengths <- this$.seqLengths;
  if (force || is.null(seqLengths)) {
    if (!onlyIfCached) {
      seqLengths <- readSeqLengths(this, ...);
      this$.seqLengths <- seqLengths;
    }
  }
  seqLengths;
})

setMethodS3("getTotalSeqLengths", "FastaReferenceFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...);
  if (is.null(seqLengths)) return(as.integer(NA));
  res <- sum(as.numeric(seqLengths));
  if (res < .Machine$integer.max) {
    res <- as.integer(res);
  }
  res;
})

setMethodS3("getSeqNames", "FastaReferenceFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...);
  if (is.null(seqLengths)) return(as.character(NA));
  names(seqLengths);
})

setMethodS3("nbrOfSeqs", "FastaReferenceFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...);
  if (is.null(seqLengths)) return(as.integer(NA));
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


###########################################################################/**
# @RdocMethod buildIndex
#
# @title "Builds an FAI index file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to @see "Rsamtools::indexFa".}
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
setMethodS3("buildIndex", "FastaReferenceFile", function(this, ..., skip=TRUE, verbose=FALSE) {
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


  verbose && enter(verbose, "Building FASTA FAI index");
  pathname <- getPathname(this);
  verbose && cat(verbose, "FASTA pathname: ", pathname);

  pathnameFAI <- sprintf("%s.fai", pathname);
  verbose && cat(verbose, "FASTA FAI pathname: ", pathnameFAI);

  pathnameFAI <- Arguments$getWritablePathname(pathnameFAI, mustNotExist=FALSE);
  if (!skip || !isFile(pathnameFAI)) {
    verbose && enter(verbose, "Building index using Rsamtools");
    require("Rsamtools") || throw("Package not loaded: Rsamtools");
    pathnameD <- indexFa(file=pathname);
    verbose && cat(verbose, "Generated file: ", pathname);
    verbose && exit(verbose);
  }
  pathnameFAI <- Arguments$getReadablePathname(pathnameFAI);

  verbose && exit(verbose);

  invisible(pathnameFAI);
}) # buildIndex()


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
#  \item{method}{A @character string specifying the algorithm to use,
#     cf. @see "bwaIndexPrefix"}
#  \item{...}{Additional arguments passed to @see "bwaIndex".}
#  \item{skip}{If @TRUE, the index files are not rebuilt if already available.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "R.filesets::GenericDataFileSet" consisting of the BWA index files.
# }
#
# \seealso{
#   Internally, @see "bwaIndex" is used.
# }
#
# @author
#*/###########################################################################
setMethodS3("buildBwaIndexSet", "FastaReferenceFile", function(this, method, ..., skip=TRUE, verbose=FALSE) {
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
  # Argument 'method':
#  method <- Arguments$getCharacter(method);
  choices <- eval(formals(bwaIndexPrefix.default)$method);
  method <- match.arg(method, choices=choices);

  # Argument 'skip':
  skip <- Arguments$getLogical(skip);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Building BWA index files");
  stopifnot(isCapableOf(aroma.seq, "bwa"));

  pathnameFA <- getPathname(this);
  verbose && cat(verbose, "FASTA reference file to be indexed: ", pathnameFA);
  verbose && cat(verbose, "Algorithm: ", method);

  # The index prefix
  prefix <- bwaIndexPrefix(pathnameFA, method=method, ...);
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

  # Read sequences information, if not already done
  n <- nbrOfSeqs(this);
  verbose && print(verbose, this);

  res <- bwaIndex(pathnameFA, indexPrefix=prefix, a=method, ..., verbose=less(verbose, 5));

  if (res != 0L) {
    throw("Failed to build BWA index. Return code: ", res);
  }

  res <- BwaIndexSet$byPrefix(prefix);

  # Sanity check
  stopifnot(!is.null(res));

  verbose && exit(verbose);

  res;
}) # buildBwaIndexSet()




###########################################################################/**
# @RdocMethod buildBowtie2IndexSet
#
# @title "Builds a Bowtie2 index files set"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to @see "bowtie2Build".}
#  \item{skip}{If @TRUE, the index files are not rebuilt if already available.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "R.filesets::GenericDataFileSet" consisting of the bowtie2 index files.
# }
#
#% \section{Benchmarking}{
#%   Examples of processing times:
#%   \itemize{
#%    \item human_g1k_v37.fasta: ~120 minutes on System A.
#%   }
#%   where:
#%   \itemize{
#%    \item 'System A': one core on a random node on a Linux 64-bit cluster.
#%    \item 'System B' is Windows 7 Pro 64-bit on a Lenovo Thinkpad X201.
#% }
#
# \seealso{
#   Internally, @see "bowtie2Build" is used.
# }
#
# @author
#*/###########################################################################
setMethodS3("buildBowtie2IndexSet", "FastaReferenceFile", function(this, ..., skip=TRUE, verbose=FALSE) {
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


  verbose && enter(verbose, "Building Bowtie2 index file");
  stopifnot(isCapableOf(aroma.seq, "bowtie2"));

  pathnameFA <- getPathname(this);
  verbose && cat(verbose, "FASTA reference file to be indexed: ", pathnameFA);
##  verbose && cat(verbose, "Algorithm: ", method);

  # The index prefix
  prefix <- bowtie2IndexPrefix(pathnameFA, ...);
  verbose && cat(verbose, "Prefix for index files: ", prefix);

  # Locate existing index files
  res <- tryCatch({
    Bowtie2IndexSet$byPrefix(prefix);
  }, error=function(ex) Bowtie2IndexSet());

  # Nothing todo?
  if (skip && isComplete(res)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && exit(verbose);
    return(res);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build bowtie2 index set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- bowtie2Build(pathnameFA, prefix, ..., verbose=less(verbose, 5));
  verbose && str(verbose, res);

  status <- attr(res, "status");
  str(status);
  if (is.null(status)) status <- 0L;
  if (status != 0L) {
    throw("Failed to build Bowtie2 index. Return code: ", status);
  }

  is <- Bowtie2IndexSet$byPrefix(prefix);

  # Sanity check
  stopifnot(!is.null(is));

  verbose && exit(verbose);

  is;
}) # buildBowtie2IndexSet()



############################################################################
# HISTORY:
# 2013-11-01
# o Now buildBowtie2IndexSet() for FastaReferenceFile supports gzip'ed
#   FASTA files.
# 2013-06-27
# o Added Rdoc comments for buildBowtie2IndexSet().
# o BUG FIX: buildBowtie2IndexSet() for FastaReferenceFile became broken
#   after updates in systemBowtie2Build().
# 2012-10-31
# o Added buildIndex() for FastaReferenceFile for building FAI index files.
# 2012-09-27
# o Added buildBowtie2IndexSet() for FastaReferenceFile.
# 2012-09-25
# o SPEEDUP: Now print() only displays sequence information, iff
#   already cached.  Otherwise, NAs are displayed.
# o Added mandatory argument 'method' to buildBwaIndexSet().
# 2012-09-24
# o Now getTotalSeqLengths() returns a numeric if the result cannot be
#   held in an integer.
# o Added buildBwaIndexSet().
# 2012-06-28
# o Created.
############################################################################
