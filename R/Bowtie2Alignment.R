###########################################################################/**
# @RdocClass Bowtie2Alignment
#
# @title "The Bowtie2Alignment class"
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
#  \item{ds}{An @see "FastqDataSet".}
#  \item{reference}{An @see "FastaReferenceFile".} of 
#  \item{tags}{Additional tags for the output data sets.}
#  \item{...}{Additional arguments passed to the aligner.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \author{Henrik Bengtsson and Pierre Neuvial}
#*/########################################################################### 
setConstructorS3("Bowtie2Alignment", function(ds=NULL, reference=NULL, tags="*", ...) {
  # Validate arguments
  if (!is.null(ds)) {
    # Argument 'ds':
    ds <- Arguments$getInstanceOf(ds, "FastqDataSet");

    # Argument 'reference':
    reference <- Arguments$getInstanceOf(ds, "FastaReferenceFile");
  } # if (!is.null(ds))

  # Arguments '...':
  args <- list(...);

  this <- extend(Object(...), "Bowtie2Alignment",
    .ds = ds,
    .reference = reference,
    .tags = tags,
    .args = args
  );

  setTags(this, tags);

  this;
})


setMethodS3("as.character", "Bowtie2Alignment", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  ds <- getInputDataSet(this);
  s <- c(s, "Input data set:");
  s <- c(s, as.character(ds));

  ref <- getReferenceFile(this);
  s <- c(s, "Aligning to reference:");
  s <- c(s, as.character(ref));
 
  class(s) <- "GenericSummary";
  s;
}, private=TRUE)


setMethodS3("getInputDataSet", "Bowtie2Alignment", function(this, ...) {
  this$.ds;
})

setMethodS3("getReferenceFile", "Bowtie2Alignment", function(this, ...) {
  this$.reference;
})

setMethodS3("getAsteriskTags", "Bowtie2Alignment", function(this, collapse=NULL, ...) {
  tags <- "bowtie2";

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, private=TRUE)


setMethodS3("getName", "Bowtie2Alignment", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
})

setMethodS3("getFlavor", "Bowtie2Alignment", function(this, ...) {
  this$.flavor;
}, protected=TRUE)


setMethodS3("getTags", "Bowtie2Alignment", function(this, collapse=NULL, ...) {
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


setMethodS3("setTags", "Bowtie2Alignment", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }
  
  this$.tags <- tags;
})

 
setMethodS3("getFullName", "Bowtie2Alignment", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "Bowtie2Alignment", function(this, ...) {
  "bamData";
})

setMethodS3("getPath", "Bowtie2Alignment", function(this, create=TRUE, ...) {
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


setMethodS3("nbrOfFiles", "Bowtie2Alignment", function(this, ...) {
  ds <- getInputDataSet(this);
  nbrOfFiles(ds);
})


setMethodS3("getOutputDataSet", "Bowtie2Alignment", function(this, ...) {
  ds <- getInputDataSet(this);
  path <- getPath(this);
  res <- BamDataSet$byPath(path, ...);
  ## TODO: Keep only those samples that exists in the input data set
  ## TODO: Assert completeness
  res;
})


setMethodS3("process", "Bowtie2Alignment", function(this, ..., force=FALSE, verbose=FALSE) {
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

  verbose && enter(verbose, "Bowtie2 alignment");
  ds <- getInputDataSet(this);
  verbose && cat(verbose, "Input data set:");
  verbose && print(verbose, ds);

  ref <- getReferenceFile(this);
  verbose && cat(verbose, "Aligning to reference:");
  verbose && print(verbose, ref);

  nbrOfFiles <- nbrOfFiles(this);
  verbose && cat(verbose, "Number of files: ", nbrOfFiles);

  outPath <- getPath(this);
  for (kk in seq(length=nbrOfFiles)) {
    df <- getFile(ds, kk);
    name <- getName(df);
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", 
                                                    kk, name, nbrOfFiles));

    # Output file
    filename <- getFilename(df);
    pathname <- Arguments$getReadablePathname(filename, path=outPath, mustExist=FALSE);

    # Nothing to do?
    if (isFile(pathname)) {
      verbose && cat(verbose, "Already aligned.");
      verbose && exit(verbose);
      next;
    }

    verbose && print(verbose, df);

    verbose && enter(verbose, "Aligning");
    dfOut <- processOne(this, df=df, reference=ref, ..., verbose=less(verbose));
    verbose && print(verbose, dfOut);
    verbose && exit(verbose);

    verbose && exit(verbose);
  } # for (kk ...)

  res <- getOutputDataSet(this, verbose=less(verbose, 1)); 

  verbose && exit(verbose);

  invisible(res);
})

setMethodS3("processOne", "Bowtie2Alignment", function(this, df, reference, ..., verbose=FALSE) {
  # Argument 'df':
  df <- Arguments$getInstanceOf(df, "FastqDataFile");

  # Argument 'reference':
  reference <- Arguments$getInstanceOf(reference, "FastaReferenceFile");


  verbose && enter(verbose, "Aligning one file");
  verbose && print(verbose, df);
  verbose && print(verbose, reference);

  fastqPathname <- getPathname(df);
  fastaPathname <- getPathname(reference);
  args <- list(pathname=fastqPathname, refPathname=fastaPathname);

  # Additional arguments
  userArgs <- list(...);
  args <- c(args, this$.args, userArgs);
  args$verbose <- less(verbose, 5);

  # Call the aligner
  bamPathname <- do.call(bowtie2, args=args);

  # Assert output
  dfOut <- BamFile(bamPathname);
  verbose && print(verbose, dfOut);
  
  verbose && exit(verbose);

  # Assert correctness of output
  dfOut <- Arguments$getInstanceOf(dfOut, "BamFile");

  dfOut;
}, private=TRUE);

############################################################################
# HISTORY:
# 2012-08-20
# o Renamed Bowtie2Aligment to Bowtie2Alignment.
# 2012-07-11
# o ROBUSTNESS: Now using do.call(bowtie2, ...) instead of 
#   do.call("bowtie2", ...) so that R CMD check can validate it
#   knowing that bowtie2() is a function.
# 2012-06-28
# o Created.
############################################################################ 
