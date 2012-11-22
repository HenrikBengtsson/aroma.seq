###########################################################################/**
# @RdocClass AbstractAlignment
#
# @title "The AbstractAlignment class"
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
#  \item{indexSet}{An @see "AbstractIndexSet".}
#  \item{tags}{Additional tags for the output data sets.}
#  \item{rgSet}{(optional) An @see "SamReadGroup" for added 
#    SAM read group to the results.}
#  \item{...}{Additional alignment arguments.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \author{Henrik Bengtsson}
#*/########################################################################### 
setConstructorS3("AbstractAlignment", function(dataSet=NULL, indexSet=NULL, tags="*", rgSet=NULL, ...) {
  # Validate arguments
  if (!is.null(dataSet)) {
    # Argument 'dataSet':
    dataSet <- Arguments$getInstanceOf(dataSet, "FastqDataSet");

    # Argument 'indexSet':
    indexSet <- Arguments$getInstanceOf(indexSet, "AbstractIndexSet");

    # Argument 'rgSet':
    if (is.null(rgSet)) {
      rgSet <- getSamReadGroup(dataSet);
    } else {
      rgSet <- Arguments$getInstanceOf(rgSet, "SamReadGroup");
    }
  } # if (!is.null(dataSet))

  # Arguments '...':
  args <- list(...);


  this <- extend(Object(), c("AbstractAlignment", uses("ParametersInterface")),
    .ds = dataSet,
    .indexSet = indexSet,
    .tags = tags,
    .rgSet = rgSet,
    .args = args
  );

  setTags(this, tags);

  this;
})


setMethodS3("as.character", "AbstractAlignment", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- sprintf("%s:", class(this)[1]);

  ds <- getInputDataSet(this);
  s <- c(s, "Input data set:");
  s <- c(s, as.character(ds));

  is <- getIndexSet(this);
  s <- c(s, "Reference index set:");
  s <- c(s, as.character(is));

  # Algorithm parameters
  s <- c(s, sprintf("Parameters: %s", getParametersAsString(this)));
  
  class(s) <- "GenericSummary";
  s;
}, protected=TRUE)


setMethodS3("getInputDataSet", "AbstractAlignment", function(this, ...) {
  this$.ds;
})

setMethodS3("getIndexSet", "AbstractAlignment", function(this, ...) {
  this$.indexSet;
})

setMethodS3("getOptionalArguments", "AbstractAlignment", function(this, ...) {
  this$.args;
}, protected=TRUE)

setMethodS3("getParameterSets", "AbstractAlignment", function(this, ...) {
  paramsList <- NextMethod("getParameterSets");
  paramsList;
}, protected=TRUE)



setMethodS3("getAsteriskTags", "AbstractAlignment", function(this, collapse=NULL, ...) {
  is <- getIndexSet(this);

  alignName <- gsub("Alignment", "", class(this)[1], fixed=TRUE);
  alignName <- tolower(alignName);
  tags <- c(alignName, getTags(is, collapse=NULL));
  tags <- unique(tags);

  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse);
  }

  tags;
}, protected=TRUE)


setMethodS3("getName", "AbstractAlignment", function(this, ...) {
  ds <- getInputDataSet(this);
  getName(ds);
})

setMethodS3("getFlavor", "AbstractAlignment", function(this, ...) {
  this$.flavor;
}, protected=TRUE)


setMethodS3("getTags", "AbstractAlignment", function(this, collapse=NULL, ...) {
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


setMethodS3("setTags", "AbstractAlignment", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags);
    tags <- trim(unlist(strsplit(tags, split=",")));
    tags <- tags[nchar(tags) > 0];
  }
  
  this$.tags <- tags;
})

 
setMethodS3("getFullName", "AbstractAlignment", function(this, ...) {
  name <- getName(this);
  tags <- getTags(this);
  fullname <- paste(c(name, tags), collapse=",");
  fullname <- gsub("[,]$", "", fullname);
  fullname;
})


setMethodS3("getRootPath", "AbstractAlignment", function(this, ...) {
  alignName <- gsub("Alignment", "", class(this)[1], fixed=TRUE);
  alignName <- tolower(alignName);
  sprintf("%sData", alignName);
}, protected=TRUE)


setMethodS3("getPath", "AbstractAlignment", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this);

  # Full name
  fullname <- getFullName(this);

  # Platform    
  ds <- getInputDataSet(this);
  platform <- "Generic";

  # The full path
  path <- filePath(rootPath, fullname, platform);

  if (create) {
    path <- Arguments$getWritablePath(path);
  } else {
    path <- Arguments$getReadablePath(path, mustExist=FALSE);
  }

  # Verify that it is not the same as the input path
  inPath <- getPath(ds);
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath);
  }

  path;
}, protected=TRUE)


setMethodS3("nbrOfFiles", "AbstractAlignment", function(this, ...) {
  ds <- getInputDataSet(this);
  length(ds);
})


setMethodS3("getOutputDataSet", "AbstractAlignment", function(this, ...) {
  ## Find all existing output data files
  path <- getPath(this);
  res <- BamDataSet$byPath(path, ...);

  ## Keep only those samples that exists in the input data set
  ds <- getInputDataSet(this);
  res <- extract(res, getFullNames(ds));
  
  ## TODO: Assert completeness
  res;
})


setMethodS3("process", "AbstractAlignment", abstract=TRUE);


############################################################################
# HISTORY:
# 2012-10-01
# o Created; extracted from BwaAlignment.R.
############################################################################ 
