###########################################################################/**
# @RdocClass BwaAlignment
#
# @title "The BwaAlignment class"
#
# \description{
#  @classhierarchy
#
#  ...
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Arguments passed to @see "AbstractAlignment".}
#  \item{indexSet}{An @see "BwaIndexSet".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \references{
#   [1] Li H. and Durbin R., \emph{Fast and accurate short read alignment
#       with Burrows-Wheeler Transform}. Bioinformatics, 2009.\cr
#   [2] Li H. and Durbin R., \emph{Fast and accurate long-read alignment
#       with Burrows-Wheeler Transform}. Bioinformatics, 2010.\cr
# }
#
# \author{Henrik Bengtsson}
#*/###########################################################################
setConstructorS3("BwaAlignment", function(..., indexSet=NULL) {
  # Validate arguments
  if (!is.null(indexSet)) {
    indexSet <- Arguments$getInstanceOf(indexSet, "BwaIndexSet");
  }

  # Arguments '...':
  args <- list(...);

  extend(AbstractAlignment(..., indexSet=indexSet), "BwaAlignment");
})


setMethodS3("getParameterSets", "BwaAlignment", function(this, ...) {
  paramsList <- NextMethod("getParameterSets");
  paramsList$aln <- getOptionalArguments(this, ...);
  paramsList;
}, protected=TRUE)



###########################################################################/**
# @RdocMethod process
#
# @title "Runs the BWA aligner"
#
# \description{
#   @get "title" on all input files.
#   The generated BAM files are all sorted and indexed.
# }
#
# @synopsis
#
# \arguments{
#  \item{...}{Not used.}
#  \item{skip}{If @TRUE, already processed files are skipped.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "BamDataSet".
# }
#
# @author
#*/###########################################################################
setMethodS3("process", "BwaAlignment", function(this, ..., skip=TRUE, force=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  asBwaParameter <- function(rg, ...) {
    if (isEmpty(rg)) {
      return(NULL);
    }

    # Validate
    if (!hasID(rg)) {
      msg <- sprintf("BWA 'samse/sampe' requires that the SAM read group has an ID. Will use default ID=1: %s", paste(as.character(rg), collapse=", "));
      warning(msg);
      rg$ID <- 1L;
    }

    rgArg <- asString(rg, fmtstr="%s:%s", collapse="\t");
    rgArg <- sprintf("@RG\t%s", rgArg);
    # Don't forget to put within quotation marks
    rgArg <- sprintf("\"%s\"", rgArg);

    rgArg <- list(r=rgArg);
  } # asBwaParameter()


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "BWA alignment");

  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  # Already done?
  todo <- findFilesTodo(this, verbose=less(verbose, 1));
  if (!force && length(todo) == 0L) {
    verbose && cat(verbose, "Already done. Skipping.");
    res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));
    verbose && exit(verbose);
    return(invisible(res));
  }

  is <- getIndexSet(this);
  verbose && cat(verbose, "Aligning using index set:");
  verbose && print(verbose, is);
  indexPrefix <- getIndexPrefix(is);

  rgSet <- this$.rgSet;
  if (!is.null(rgSet)) {
    verbose && cat(verbose, "Assigning SAM read group:");
    verbose && print(verbose, rgSet);
    validate(rgSet);
  }

  paramsList <- getParameterSets(this);
  verbose && printf(verbose, "Additional BWA arguments: %s\n", getParametersAsString(this));

  verbose && cat(verbose, "Number of files: ", length(ds));

  isPaired <- isPaired(this);
  if (isPaired) {
    pairs <- getFilePairs(ds);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Apply aligner to each of the FASTQ files
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  dsApply(ds[todo], FUN=function(df, paired=FALSE, indexPrefix, rgSet, paramsList, path, ...., skip=TRUE, verbose=FALSE) {
    R.utils::use("R.utils, aroma.seq");

    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Validate arguments
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    # Argument 'df':
    df <- Arguments$getInstanceOf(df, "FastqDataFile");

    # Argument 'paired':
    paired <- Arguments$getLogical(paired);

    # Argument 'indexPrefix':
    indexPrefix <- Arguments$getCharacter(indexPrefix);

    # Argument 'skip':
    skip <- Arguments$getLogical(skip);

    # Argument 'path':
    path <- Arguments$getWritablePath(path);

    # Argument 'verbose':
    verbose <- Arguments$getVerbose(verbose);
    if (verbose) {
      pushState(verbose);
      on.exit(popState(verbose));
    }


    verbose && enter(verbose, "BWA alignment of one sample");

    if (isPaired) {
      pathnameFQ <- sapply(list(R1=df, R2=getMateFile(df)), FUN=getPathname);
      verbose && cat(verbose, "FASTQ R1 pathname: ", pathnameFQ[1L]);
      verbose && cat(verbose, "FASTQ R2 pathname: ", pathnameFQ[2L]);
      stop("Paired-end BWA alignment is not yet supported.");
    } else {
      pathnameFQ <- getPathname(df);
      verbose && cat(verbose, "FASTQ pathname: ", pathnameFQ);
    }

    # The SAI and SAM files to be generated
    fullname <- getFullName(df);
    filename <- sprintf("%s.sai", fullname);
    pathnameSAI <- Arguments$getWritablePathname(filename, path=path);
    filename <- sprintf("%s.sam", fullname);
    pathnameSAM <- Arguments$getWritablePathname(filename, path=path);
    verbose && cat(verbose, "SAM pathname: ", pathnameSAM);
    filename <- sprintf("%s.bam", fullname);
    pathnameBAM <- Arguments$getWritablePathname(filename, path=path);
    verbose && cat(verbose, "BAM pathname: ", pathnameBAM);

    # Nothing to do?
    done <- (skip && isFile(pathnameBAM));
    if (done) {
      verbose && cat(verbose, "Already aligned. Skipping");
    } else {
      verbose && print(verbose, df);

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (a) Generate SAI file via BWA aln
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!isFile(pathnameSAI)) {
        args <- list(pathnameFQ, indexPrefix=indexPrefix,
                     pathnameD=pathnameSAI);
        args <- c(args, paramsList$aln);
        args$verbose <- less(verbose, 5);
        res <- do.call(bwaAln, args=args);
        verbose && cat(verbose, "System result code: ", res);
      }

      # Sanity check
      stopifnot(isFile(pathnameSAI));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (b) Generate SAM file via BWA samse
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!isFile(pathnameSAM)) {
        # Extract sample-specific read group
        rgII <- getSamReadGroup(df);
        if (length(rgSet) > 0L) {
          rgII <- merge(rgSet, rgII);
        }
        verbose && cat(verbose, "Writing SAM Read Groups:");
        verbose && print(verbose, rgII);
        verbose && cat(verbose, "BWA 'samse' parameter:");
        rgArg <- asBwaParameter(rgII);
        verbose && print(verbose, rgArg);

        args <- list(pathnameSAI=pathnameSAI, pathnameFQ=pathnameFQ,
                     indexPrefix=indexPrefix, pathnameD=pathnameSAM);
        args <- c(args, rgArg);
        args$verbose <- less(verbose, 5);
        res <- do.call(bwaSamse, args=args);
        verbose && cat(verbose, "System result code: ", res);

        # BWA 'samse' can generate empty SAM files
        if (isFile(pathnameSAM)) {
          if (file.info(pathnameSAM)$size == 0L) {
            verbose && cat(verbose, "Removing empty SAM file falsely created by BWA: ", pathnameSAM);
            file.remove(pathnameSAM);
          }
        }
      }
      # Sanity check
      stopifnot(isFile(pathnameSAM));

      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      # (c) Generates a (sorted and indexed) BAM file from SAM file
      # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      if (!isFile(pathnameBAM)) {
        sf <- SamDataFile(pathnameSAM);
        bf <- convertToBam(sf, verbose=less(verbose, 5));
        verbose && print(verbose, pathnameBAM);
      }
      # Sanity check
      stopifnot(isFile(pathnameBAM));
    } # if (done)

    verbose && exit(verbose);

    invisible(list(pathnameFQ=pathnameFQ, pathnameSAM=pathnameSAM, pathnameBAM=pathnameBAM));
  }, paired=isPaired(this), indexPrefix=indexPrefix, rgSet=rgSet, paramsList=paramsList, path=getPath(this), skip=skip, verbose=verbose) # dsApply()

  res <- getOutputDataSet(this, onMissing="error", verbose=less(verbose, 1));

  verbose && exit(verbose);

  invisible(res);
})



############################################################################
# HISTORY:
# 2013-11-11
# o SPEEDUP: Now Bowtie2Alignment and BwaAlignment skips already processed
#   items much faster and if all are done, even quicker.
# 2013-08-31
# o Now process() for BwaAlignment utilizes dsApply().
# 2012-10-21
# o Added argument 'drop' to getParameters().
# o Now a default 'ID' tag is added to the SAM read group if missing.
# 2012-10-01
# o Now process() BwaAlignment write SAM read groups, iff given.
# o Now BwaAlignment inherits from AbstractAlignment.
# 2012-09-28
# o Added support for argument 'rgSet' to BwaAlignment().
# 2012-09-25
# o Now process() passes additional arguments to bwaAln().
# o Created from Bowtie2Alignment.R.
############################################################################
