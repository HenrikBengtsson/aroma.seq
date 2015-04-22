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
# \section{Compression}{
#  Currently, the package only supports non-compressed FASTA files.
# }
#
# \section{Filenames}{
#  Currently, FASTA files with commas in their filenames should be avoided
#  because they are not supported by Bowtie2.
# }
#
# \section{On different FASTA files}{
#  Some FASTA files represent nucleotides as
#  upper-case letters A, C, G, T (and N), whereas others as
#  lower-case letters a, c, g, t (and n).
#  The format of the sequence names differ between FASTA file
#  (even for the same genome/organism), e.g. 'chr1' versus '1'.
# }
#
# @author "HB"
#
# \seealso{
#   [1] \url{http://www.wikipedia.org/wiki/FASTA_format}
# }
#*/###########################################################################
setConstructorS3("FastaReferenceFile", function(...) {
  extend(GenericDataFile(...), "FastaReferenceFile",
    .seqLengths=NULL
  );
})

setMethodS3("as.character", "FastaReferenceFile", function(x, ...) {
  this <- x
  s <- NextMethod("as.character")
  n <- nbrOfSeqs(this)
  s <- c(s, sprintf("Total sequence length: %s", pi3(getTotalSeqLengths(this))))
  s <- c(s, sprintf("Number of sequences: %d", n))
  s <- c(s, sprintf("Sequence names: [%d] %s", n, hpaste(getSeqNames(this))))
  s
}, protected=TRUE)


setMethodS3("getDefaultFullName", "FastaReferenceFile", function(this, ...) {
  name <- NextMethod("getDefaultFullName");
  name <- gsub("[.](fasta|fa)$", "", name, ignore.case=TRUE);
  name;
}, protected=TRUE)


setMethodS3("getOrganism", "FastaReferenceFile", function(this, ...) {
  path <- getPath(this);
  organism <- basename(path);
  organism;
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
  if (is.null(seqLengths)) return(NA_integer_);
  res <- sum(as.numeric(seqLengths));
  if (res < .Machine$integer.max) {
    res <- as.integer(res);
  }
  res;
})

setMethodS3("getSeqNames", "FastaReferenceFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...);
  if (is.null(seqLengths)) return(NA_character_);
  names <- names(seqLengths);

  # From http://en.wikipedia.org/wiki/FASTA_format:
  # "The **word** following the ">" symbol is the identifier of the
  # sequence, and the rest of the line is the description [...]", e.g.
  # >I dna:chromosome chromosome:EF4:I:1:230218:1 REF
  # => ID/name: 'I'
  # => Description: 'dna:chromosome chromosome:EF4:I:1:230218:1 REF'
  names <- gsub(" .*", "", names);

  names;
})

setMethodS3("nbrOfSeqs", "FastaReferenceFile", function(this, ...) {
  seqLengths <- getSeqLengths(this, ...);
  if (is.null(seqLengths)) return(NA_integer_);
  length(seqLengths);
})


# \seealso{
#   Internally, \code{fasta.info()} of \pkg{Biostrings} is used.
# }
setMethodS3("readSeqLengths", "FastaReferenceFile", function(this, force=FALSE, ...) {
  pathname <- getPathname(this);

  # Check for cached results
  dirs <- c("aroma.seq", getOrganism(this));
  key <- list(method="readSeqLengths", class=class(this), pathname=pathname);
  seqLengths <- loadCache(key=key, dirs=dirs);
  if (!force && !is.null(seqLengths)) {
    return(seqLengths);
  }

  # Read FASTA file
  seqLengths <- Biostrings::fasta.info(pathname);

  # Cache
  saveCache(seqLengths, key=key, dirs=dirs);

  seqLengths;
}, private=TRUE)




###########################################################################/**
# @RdocMethod byOrganism
# @aliasmethod findByOrganism
#
# @title "Locates a FASTA file by organism"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{organism}{A @character string specifying for which organism a
#    file should be retrieved.}
#  \item{tags}{(not used) A @character @vector.}
#  \item{prefix}{(optional) A @character string specifying an optional
#    regular expression prefix to be prepended to \code{pattern} when
#    searching for the file.}
#  \item{pattern}{A @character string specifying a regular expression for
#    the file to be located.}
#  \item{...}{Additional arguments passed to the constructor of
#    @see "FastaReferenceFile" when instantiating the object.}
# }
#
# \value{
#   Returns a @see "FastaReferenceFile".
# }
#
# \seealso{
#   @seeclass
# }
#
# @author
#*/###########################################################################
setMethodS3("findByOrganism", "FastaReferenceFile", function(static, organism, tags=NULL, prefix=NULL, pattern="[.](fa|fasta)(|[.]gz)$", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'organism':
  organism <- Arguments$getOrganism(organism);

  # Argument 'prefix':
  if (!is.null(prefix)) {
    prefix <- Arguments$getRegularExpression(prefix);
  }


  args <- list(pattern=pattern);

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Search in annotationData/organisms/<organism>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create the fullname
  fullname <- paste(c(organism, tags), collapse=",");

  # Extract the name and the tags
  parts <- unlist(strsplit(fullname, split=",", fixed=TRUE));
  organism <- parts[1L];
  tags <- parts[-1L];

  # Search for "organisms/<organism>/<prefix>.*[.](fa|fasta)$" files
  patternS <- pattern;
  if (!is.null(prefix)) patternS <- sprintf("%s.*%s", prefix, patternS);
  args <- list(
    set="organisms",
    name=organism,
    pattern=patternS,
    ...
  );
  pathname <- do.call("findAnnotationData", args=args);

  # If not found, look for Windows shortcuts
  if (is.null(pathname)) {
    # Search for a Windows shortcut
    args$pattern <- sprintf("%s[.]lnk$", args$pattern)
    pathname <- do.call("findAnnotationData", args=args);
    if (!is.null(pathname)) {
      # ..and expand it
      pathname <- Arguments$getReadablePathname(pathname, mustExist=FALSE);
      if (!isFile(pathname))
        pathname <- NULL;
    }
  }

  pathname;
}, static=TRUE, protected=TRUE) # findByOrganism()


setMethodS3("byOrganism", "FastaReferenceFile", function(static, organism, ...) {
  # Argument 'organism':
  organism <- Arguments$getOrganism(organism);

  # Locate FASTA file
  pathname <- findByOrganism(static, organism, ...);
  if (length(pathname) == 0L)
    throw("Failed to located FASTA reference file for organism: ", organism);

  # Allocate object
  res <- newInstance(static, pathname, ..., .onUnknownArgs="ignore");

  # Validate
  organismR <- getOrganism(res);
  if (organismR != organism) {
    throw(sprintf("The located %s (%s) specifies an organism different from the requested one: %s != %s", class(res)[1L], getPathname(res), sQuote(organismR), sQuote(organism)));
  }
  res;
}, static=TRUE) # byOrganism()





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
#   Returns the pathname to the FASTA index file.
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
    use("Rsamtools")
    pathnameD <- indexFa(file=pathname);
    verbose && cat(verbose, "Generated file: ", pathname);
    verbose && exit(verbose);
  }
  pathnameFAI <- Arguments$getReadablePathname(pathnameFAI);

  verbose && exit(verbose);

  invisible(pathnameFAI);
}) # buildIndex()




###########################################################################/**
# @RdocMethod buildDictionary
#
# @title "Builds a DICT sequence dictionary file"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Additional arguments passed to @see "systemPicard".}
#  \item{skip}{If @TRUE, the dictionary is not rebuilt if already available.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns the pathname to the DICT file.
# }
#
# \references{
#  [1] Geraldine van der Auwera,
#      How can I prepare a FASTA file to use as reference?,
#      GATK Forum, Sept 2013.
#      \url{http://gatkforums.broadinstitute.org/discussion/1601/how-can-i-prepare-a-fasta-file-to-use-as-reference}
# }
#
# @author
#*/###########################################################################
setMethodS3("buildDictionary", "FastaReferenceFile", function(this, ..., skip=TRUE, verbose=FALSE) {
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


  verbose && enter(verbose, "Building FASTA sequence dictionary");
  pathname <- getPathname(this);
  verbose && cat(verbose, "FASTA pathname: ", pathname);

  path <- getPath(this);
  filename <- sprintf("%s.dict", getFullName(this));
  pathnameDICT <- file.path(path, filename);
  verbose && cat(verbose, "FASTA DICT pathname: ", pathnameDICT);

  pathnameDICT <- Arguments$getWritablePathname(pathnameDICT, mustNotExist=FALSE);
  if (!skip || !isFile(pathnameDICT)) {
    verbose && enter(verbose, "Building sequence dictionary using Picard");

    # Emulate overwriting
    if (isFile(pathnameDICT)) file.remove(pathnameDICT);

    # Write to temporary file
    pathnameDICTT <- pushTemporaryFile(pathnameDICT);

    res <- systemPicard(command="CreateSequenceDictionary", R=pathname, O=pathnameDICTT, verbose=less(verbose, 50));
    verbose && cat(verbose, "System result:");
    verbose && str(verbose, res);

    pathnameDICT <- popTemporaryFile(pathnameDICTT);
    verbose && cat(verbose, "Generated file: ", pathnameDICT);
    verbose && exit(verbose);
  }
  pathnameDICT <- Arguments$getReadablePathname(pathnameDICT);

  verbose && exit(verbose);

  invisible(pathnameDICT);
}) # buildDictionary()




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
#  \item{method}{A @character string specifying the algorithm to use for
#     building the index set. All methods gives identical results [1].
#     The default is such that it can handle also large genomes, including
#     the human genome.}
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
# \references{
#  [1] Thread \emph{bwa index option bwtsw}, SEQanswers, 2010-07-13.
#      \url{http://seqanswers.com/forums/showthread.php?t=5921}\cr
#  [2] Edwards Bioinformatics Lab,
#      \emph{How to create a database for BWA and BWA-SW}, 2013.
#      \url{http://edwards.sdsu.edu/research/index.php/robert/282-how-to-create-a-database-for-bwa-and-bwa-sw} \cr
#  [3] Henrik Bengtsson, \emph{bwa index -a is: Details on database 2GB limit?},
#      bwa-help thread on 2013-11-18.
#      \url{https://sourceforge.net/mailarchive/message.php?msg_id=31649355}\cr
# }
#
# @author
#*/###########################################################################
setMethodS3("buildBwaIndexSet", "FastaReferenceFile", function(this, method=c("bwtsw", "is"), ..., skip=TRUE, verbose=FALSE) {
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
  method <- match.arg(method);

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

  # Sanity check
  if (method == "is") {
    verbose && enter(verbose, "Asserting that index can be build with method 'is'");
    size <- getFileSize(this);
    verbose && printf(verbose, "FASTA filesize: %.f bytes\n", size);
    maxNbrOfBases <- 2e9;
    if (size > maxNbrOfBases) {
      nbrOfBases <- getTotalSeqLengths(this);
      verbose && printf(verbose, "Number of bases in FASTA reference: %.f\n", nbrOfBases);
      if (nbrOfBases > maxNbrOfBases) {
        throw(sprintf("Cannot build BWA index with method 'is' (consider using 'bwtsw' instead).  There are too many bases in FASTA file (%s): %.0f > %.0f", sQuote(pathnameFA), nbrOfBases, maxNbrOfBases));
      }
    }
    verbose && exit(verbose);
  }

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

  # Locate existing index set
  is <- tryCatch({
    Bowtie2IndexSet$byPrefix(prefix);
  }, error=function(ex) Bowtie2IndexSet());

  # Nothing todo?
  if (skip && isComplete(is)) {
    verbose && cat(verbose, "Already done. Skipping.");
    verbose && print(verbose, is)
    # Assert compatibility
    isCompatibleWith(this, is, mustWork=TRUE, verbose=less(verbose, 50))
    verbose && exit(verbose);
    return(is);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Build bowtie2 index set
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  res <- bowtie2Build(pathnameFA, prefix, ..., verbose=less(verbose, 5));
  status <- attr(res, "status"); if (is.null(status)) status <- 0L;
  verbose && cat(verbose, "Results:");
  verbose && str(verbose, res);
  verbose && cat(verbose, "Status:");
  verbose && str(verbose, status);
  if (status != 0L) {
    throw("Failed to build Bowtie2 index. Return code: ", status);
  }

  is <- Bowtie2IndexSet$byPrefix(prefix);

  # Sanity check
  stopifnot(!is.null(is));

  # Assert compatibility
  isCompatibleWith(this, is, mustWork=TRUE, verbose=less(verbose, 50))

  verbose && exit(verbose);

  is;
}) # buildBowtie2IndexSet()


############################################################################
# HISTORY:
# 2014-08-23
# o ROBUSTNESS: Now buildBowtie2IndexSet() asserts that the returned
#   index set is compatible with the FASTA file.
# 2014-08-11
# o BUG FIX: getSeqNames() for FastaReferenceFile would return the
#   sequence description in addition to the ID as part of the name.
# 2014-04-16
# o SPEEDUP: Now readSeqLengths() for FastaReferenceFile memoizes results.
# 2014-04-13
# o Added buildDictionary() for FastaReferenceFile.
# 2014-02-25
# o Now static byOrganism() no longer passes '...'.
# 2014-02-20
# o Now findByOrganism() for FastaReferenceFile also locates gzipped files.
# o Analogously to FastqDataFile, added getDefaultFullName() for
#   FastaReferenceFile so <fullname>.fasta.gz is properly handled.
#   Should ideally handled by R.filesets.
# 2014-01-25
# o Now static byOrganism() passes '...' also to the constructor.
# o DOCUMENTATION: Added help for byOrganism().
# 2013-11-19
# o ROBUSTNESS: Now FastaReferenceFile$byOrganism() asserts that the
#   returned FASTA file specifies the requested organism.
# 2013-11-17
# o Now buildBwaIndexSet(..., method="is") checks for maximum size of
#   reference genome (2GB) and gives an informative error if the FASTA
#   file is greater.
# o Added argument 'prefix' to findByOrganism() for FastaReferenceFile.
# 2013-11-10
# o Added getOrganism().
# 2013-11-09
# o Added FastaReferenceFile$byOrganism().
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
