###########################################################################/**
# @set "class=Bowtie2IndexSet"
# @RdocMethod buildTopHat2TranscriptomeIndexSet
# @alias buildTopHat2TranscriptomeIndexSet
#
# @title "Calls TopHat to build a transcriptome index"
#
# \description{
#  @get "title".
# }
#
# \arguments{
#   \item{gtf}{GtfDataFile to be indexed}
#   \item{outPath}{(optional) Output directory for index and log file.}
#   \item{...}{Arguments passed to tophat().}
#   \item{skip}{If @TRUE, the index files are not rebuilt if already available.}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# \references{
#  [1] TopHat - A spliced read mapper for RNA-Seq, 2015.
#      \url{http://ccb.jhu.edu/software/tophat/}
# }
#
# @keyword internal
#*/###########################################################################
setMethodS3("buildTopHat2TranscriptomeIndexSet", "Bowtie2IndexSet", function(this, gtf, fa=NULL, outPath=NULL, ..., skip=TRUE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  hasCommas <- function(pathnames, ...) {
    (regexpr(",", pathnames, fixed=TRUE) != -1L);
  } # hasCommas()

  assertNoCommas <- function(pathnames, ...) {
    stopifnot(!any(hasCommas(pathnames)));
  } # assertNoCommas()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'gtf':
  pathnameGTF <- getPathname(gtf)
  pathnameGTF <- Arguments$getReadablePathname(pathnameGTF)
  assertNoCommas(getFullName(gtf))

  # Argument 'fa':
  if (!is.null(fa)) {
    fa <- Arguments$getInstanceOf(fa, "FastaReferenceFile");
  }

  # Argument 'outPath':
  if (is.null(outPath)) outPath <- file.path(getParent(getPath(this)), "tophat2")
  outPath <- Arguments$getReadablePath(outPath, mustExist=FALSE)
  if (!isDirectory(outPath)) {
    outPath <- Arguments$getWritablePath(outPath)
  }

  # Argument 'skip':
  skip <- Arguments$getLogical(skip)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose), add=TRUE)
  }


  verbose && enter(verbose, "Building TopHat transcriptome index")
  verbose && print(verbose, "Bowtie2 index set:")
  verbose && print(verbose, this)

  verbose && cat(verbose, "GTF:")
  verbose && print(verbose, gtf)

  # The (path) prefix used for the transcripome index set
  prefixName <- getFullName(gtf)
  prefix <- file.path(outPath, prefixName)
  verbose && cat(verbose, "TopHat2 transcriptome index set prefix: ", prefix)

  # Check for an existing transcriptome and Bowtie2 index set
  # cf. http://ccb.jhu.edu/software/tophat/manual.shtml#tbuild
  tis <- tryCatch({
    Bowtie2IndexSet$byPrefix(prefix)
  }, error=function(ex) Bowtie2IndexSet())

  # Nothing todo?
  if (skip && isComplete(tis)) {
    verbose && cat(verbose, "Transcriptome indexing already done. Skipping.")
    verbose && print(verbose, tis)

    # Assert compatibility
    isCompatibleWith(tis, gtf, mustWork=TRUE, verbose=less(verbose, 50))

    verbose && exit(verbose)

    return(tis)
  }

  # Assert compatibility
  isCompatibleWith(this, gtf, mustWork=TRUE, verbose=less(verbose, 50))


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Call TopHat executable
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Assert TopHat and that it is of a sufficient version
  stopifnot(isCapableOf(aroma.seq, "tophat2"))
  verT <- attr(findTopHat2(), "version")
  if (verT < '2.0.10') throw("TopHat version >= 2.0.10 required")

  # Link to existing FASTA file (avoids having tophat2 to recreate it)
  if (is.null(fa)) fa <- getFastaReferenceFile(this);
  verbose && cat(verbose, "FASTA:");
  verbose && print(verbose, fa);
  linkTo(fa, path=getPath(this));

  # Output to temporary directory
  outPathT <- sprintf("%s.tmp", outPath);

  # (Pre-existing) index for the reference genome
  res <- tophat(getIndexPrefix(this), gtf=pathnameGTF, outPath=outPathT,
                optionsVec=c("--transcriptome-index"="."),
                ..., verbose=less(verbose, 10))
  status <- attr(res, "status"); if (is.null(status)) status <- 0L;
  verbose && cat(verbose, "Results:");
  verbose && str(verbose, res);
  verbose && cat(verbose, "Status: ", status);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Cleanup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove local tmp/ directory, if exists and empty
  pathT <- file.path(outPathT, "tmp");
  if (isDirectory(pathT) && length(listDirectory(pathT, private=TRUE)) == 0L) {
    removeDirectory(pathT);
  }

  # Rename any remaining local directories
  dirs <- listDirectory(outPathT, private=TRUE, fullNames=TRUE);
  dirs <- dirs[isDirectory(dirs)];
  for (pathT in dirs) {
    dirTT <- sprintf("%s.%s", getFullName(gtf), basename(pathT));
    pathTT <- file.path(outPathT, dirTT);
    file.rename(pathT, pathTT);
  }

  # Try to setup index set
  prefixT <- file.path(outPathT, prefixName)
  tis <- Bowtie2IndexSet$byPrefix(prefixT)
  # Assert compatibility
  isCompatibleWith(tis, gtf, mustWork=TRUE, verbose=less(verbose, 50))

  # Move all related files and directories one by one to final destination
  pattern <- sprintf("^%s[.]", prefixName);
  for (filename in listDirectory(path=outPathT, pattern=pattern)) {
    pathnameS <- file.path(outPathT, filename);
    pathnameD <- file.path(outPath, filename);
    renameFile(pathnameS, pathnameD, overwrite=FALSE);
  }

  # Assert that all files have been moved
  stopifnot(length(listDirectory(outPathT, pattern=pattern)) == 0L);

  # Remove temporary directory, if empty
  if (length(listDirectory(outPathT, private=TRUE)) == 0L) {
    removeDirectory(outPathT);
  }



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Locate index set to return
  tis <- Bowtie2IndexSet$byPrefix(prefix)
  verbose && print(verbose, tis)

  # Assert compatibility
  isCompatibleWith(tis, gtf, mustWork=TRUE, verbose=less(verbose, 50))

  verbose && exit(verbose)

  tis
}) # buildTopHat2TranscriptomeIndexSet()
