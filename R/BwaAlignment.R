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


setMethodS3("getParameters", "BwaAlignment", function(this, which=c("aln", "samse"), ..., drop=TRUE) {
  params <- list();
  params$aln <- getOptionalArguments(this, ...);

  keep <- na.omit(match(which, names(params)));
  params <- params[which];

  if (length(params) == 1L) {
    params <- params[[1L]];
  }

  params;
})

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
      msg <- sprintf("BWA 'samse/sampe' requires that the SAM read group has an ID. Will use default ID=1: %s", rgArg);
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

  
  paramsList <- getParameters(this);
  if (verbose) {
    for (key in names(paramsList)) {
      verbose && printf(verbose, "Additional BWA '%s' arguments:\n", key);
      verbose && str(verbose, paramsList[[key]]);
    }
  }


  nbrOfFiles <- nbrOfFiles(this);
  verbose && cat(verbose, "Number of files: ", nbrOfFiles);

  outPath <- getPath(this);
  for (kk in seq(length=nbrOfFiles)) {
    df <- getFile(ds, kk);
    name <- getName(df);
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", 
                                                    kk, name, nbrOfFiles));

    pathnameFQ <- getPathname(df);
    verbose && cat(verbose, "FASTQ pathname: ", pathnameFQ);

    # The SAI and SAM files to be generated
    fullname <- getFullName(df);
    filename <- sprintf("%s.sai", fullname);
    pathnameSAI <- Arguments$getWritablePathname(filename, path=outPath);
    filename <- sprintf("%s.sam", fullname);
    pathnameSAM <- Arguments$getWritablePathname(filename, path=outPath);
    verbose && cat(verbose, "SAM pathname: ", pathnameSAM);
    filename <- sprintf("%s.bam", fullname);
    pathnameBAM <- Arguments$getWritablePathname(filename, path=outPath);
    verbose && cat(verbose, "BAM pathname: ", pathnameBAM);

    # Nothing to do?
    if (skip && isFile(pathnameBAM)) {
      verbose && cat(verbose, "Already aligned. Skipping");
      verbose && exit(verbose);
      next;
    }

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
    # (c) Generate BAM file from SAM file
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    if (!isFile(pathnameBAM)) {
      sf <- SamDataFile(pathnameSAM);
      bf <- convertToBamDataFile(sf, verbose=less(verbose, 5));
      print(pathnameBAM);
    }
    # Sanity check
    stopifnot(isFile(pathnameBAM));

    verbose && exit(verbose);
  } # for (kk ...)

  res <- getOutputDataSet(this, verbose=less(verbose, 1)); 

  verbose && exit(verbose);

  invisible(res);
})



############################################################################
# HISTORY:
# 2012-10-01
# o Now process() BwaAlignment write SAM read groups, iff given.
# o Now BwaAlignment inherits from AbstractAlignment.
# 2012-09-28
# o Added support for argument 'rgSet' to BwaAlignment().
# 2012-09-25
# o Now process() passes additional arguments to bwaAln().
# o Created from Bowtie2Alignment.R.
############################################################################ 
