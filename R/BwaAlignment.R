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
#  \item{dataSet}{An @see "FastqDataSet".}
#  \item{indexSet}{An @see "BwaIndexSet".}
#  \item{tags}{Additional tags for the output data sets.}
#  \item{rgSet}{(optional) An @see "SamReadGroup" for added 
#    SAM read group to the results.}
#  \item{...}{Additional BWA 'aln' arguments.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \author{Henrik Bengtsson and Pierre Neuvial}
#*/########################################################################### 
setConstructorS3("BwaAlignment", function(dataSet=NULL, indexSet=NULL, tags="*", rgSet=NULL, ...) {
  # Validate arguments
  if (!is.null(dataSet)) {
    # Argument 'dataSet':
    dataSet <- Arguments$getInstanceOf(dataSet, "FastqDataSet");

    # Argument 'indexSet':
    if (inherits(indexSet, "FastaReferenceFile")) {
      throw("Argument 'indexSet' should be a BwaIndexSet object (not ", class(indexSet)[1], "). Use buildBwaIndexSet().");
    }
    indexSet <- Arguments$getInstanceOf(indexSet, "BwaIndexSet");

    # Argument 'rgSet':
    if (!is.null(rgSet)) {
      rgSet <- Arguments$getInstanceOf(rgSet, "SamReadGroup");
    }
  } # if (!is.null(dataSet))

  # Arguments '...':
  args <- list(...);

  this <- extend(Object(), "BwaAlignment",
    .ds = dataSet,
    .indexSet = indexSet,
    .tags = tags,
    .rgSet = rgSet,
    .args = args
  );

  setTags(this, tags);

  this;
})


setMethodS3("as.character", "BwaAlignment", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  ds <- getInputDataSet(this);
  s <- c(s, "Input data set:");
  s <- c(s, as.character(ds));

  is <- getIndexSet(this);
  s <- c(s, "Reference index set:");
  s <- c(s, as.character(is));

  # Additional arguments
  paramsList <- getParametersAsString(this, collapse=", "); 
  for (key in names(paramsList)) {
    params <- paramsList[[key]];
    s <- c(s, sprintf("Additional '%s' parameters: %s", key, params));
  }
  
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getInputDataSet", "BwaAlignment", function(this, ...) {
  this$.ds;
})

setMethodS3("getIndexSet", "BwaAlignment", function(this, ...) {
  this$.indexSet;
})

setMethodS3("getOptionalArguments", "BwaAlignment", function(this, ...) {
  this$.args;
})

setMethodS3("getParameters", "BwaAlignment", function(this, which=c("aln", "samse"), ..., drop=TRUE) {
  params <- list();
  params$aln <- getOptionalArguments(this, ...);

  rgSet <- this$.rgSet;
  if (!is.null(rgSet)) {
    rgArg <- asString(rgSet, fmtstr="%s:%s", collapse="\t");
    rgArg <- sprintf("@RG\t%s", rgArg);
    # Don't forget to put within quotation marks
    rgArg <- sprintf("\"%s\"", rgArg);
    rgArg <- list(r=rgArg);
    params$samse <- rgArg;
  }

  params <- params[which];

  if (length(params) == 1L) {
    params <- params[[1L]];
  }

  params;
})


setMethodS3("getParametersAsString", "BwaAlignment", function(this, ..., collapse=NULL, drop=TRUE) {
  paramsList <- getParameters(this, drop=FALSE, ...);
  paramsList <- lapply(paramsList, FUN=function(params) {
    params <- trim(capture.output(str(params)))[-1];
    params <- gsub("^[$][ ]*", "", params);
    params <- gsub(" [ ]*", " ", params);
    params <- gsub("[ ]*:", ":", params);
    if (!is.null(collapse)) {
      params <- paste(params, collapse=collapse);
    }
  });

  if (length(paramsList) == 1L) {
    paramsList <- paramsList[[1L]];
  }

  paramsList;
})


setMethodS3("getAsteriskTags", "BwaAlignment", function(this, collapse=NULL, ...) {
  is <- getIndexSet(this);
  tags <- c("bwa", getTags(is, collapse=NULL));
  tags <- unique(tags);

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, private=TRUE)


setMethodS3("getName", "BwaAlignment", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
})

setMethodS3("getFlavor", "BwaAlignment", function(this, ...) {
  this$.flavor;
}, protected=TRUE)


setMethodS3("getTags", "BwaAlignment", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input data set
  ds <- getInputDataSet(this);
  tags <- getTags(ds, collapse=collapse);

  # Get class-specific tags
  tags <- c(tags, this$.tags);

  # Update default tags
  tags[tags == "*"] <- getAsteriskTags(this, collapse=",");

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  } else {
    tags <- unlist(strsplit(tags, split=","));
  }

  if (length(tags) == 0) {
    tags <- NULL;
  }

  tags;
})


setMethodS3("setTags", "BwaAlignment", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }
  
  this$.tags <- tags;
})

 
setMethodS3("getFullName", "BwaAlignment", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "BwaAlignment", function(this, ...) {
  "bwaData";
})

setMethodS3("getPath", "BwaAlignment", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Platform    
  ds <- getInputDataSet(this);
  platform <- "Generic";

  # The full path
  path <- filePath(rootPath, fullname, platform, expandLinks="any");

  # Verify that it is not the same as the input path
  ds <- getInputDataSet(this);
  inPath <- getPath(ds);
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  # Create path?
  if (create) {
    if (!isDirectory(path)) {
      mkdirs(path);
      if (!isDirectory(path))
        throw("Failed to create output directory: ", path);
    }
  }

  path;
})


setMethodS3("nbrOfFiles", "BwaAlignment", function(this, ...) {
  ds <- getInputDataSet(this);
  nbrOfFiles(ds);
})


setMethodS3("getOutputDataSet", "BwaAlignment", function(this, ...) {
  ## Find all existing output data files
  path <- getPath(this);
  res <- BamDataSet$byPath(path, ...);

  ## Keep only those samples that exists in the input data set
  ds <- getInputDataSet(this);
  res <- extract(res, getFullNames(ds));
  
  ## TODO: Assert completeness
  res;
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
  # 
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  asBwaParameter <- function(rg, ...) {
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


  units <- NULL;

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
    # Extract sample-specific read group
    # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    rgII <- extractSamReadGroup(df);
    rgT <- merge(rgSet, rgII);
    verbose && cat(verbose, "Writing SAM Read Groups:");
    verbose && print(verbose, rgT);
    verbose && cat(verbose, "BWA 'samse' parameter:");
    rgArg <- asBwaParameter(rgT);
    verbose && print(verbose, rgArg);


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
      args <- list(pathnameSAI=pathnameSAI, pathnameFQ=pathnameFQ, 
                   indexPrefix=indexPrefix, pathnameD=pathnameSAM);
      args <- c(args, rgArg);
      args$verbose <- less(verbose, 5);
      res <- do.call(bwaSamse, args=args);
      verbose && cat(verbose, "System result code: ", res);
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



setMethodS3("asBwaString", "SamReadGroup", function(this, ...) {
  res <- asString(this, fmtstr="%s:%s", collapse="\t", ...)
  res <- sprintf("@RG\t%s", res);
  res;
})


############################################################################
# HISTORY:
# 2012-09-28
# o Added support for argument 'rgSet' to BwaAlignment().
# 2012-09-25
# o Now process() passes additional arguments to bwaAln().
# o Created from Bowtie2Alignment.R.
############################################################################ 
