###########################################################################/**
# @RdocClass AromaSeqTransform
#
# @title "The AromaSeqTransform class"
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
#  \item{dataSet}{An @see "GenericDataFileSet".}
#  \item{tags}{Tags appended to the output data sets.}
#  \item{flavor}{An optional @character string.}
#  \item{...}{Additional named parameters to be recorded.}
#  \item{.className}{A @character string specifying what class of
#   data sets to accept.}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setConstructorS3("AromaSeqTransform", function(dataSet=NULL, tags="*", flavor=NULL, ..., .className="GenericDataFileSet") {
  # Validate arguments
  if (!is.null(dataSet)) {
    # Argument 'dataSet':
    dataSet <- Arguments$getInstanceOf(dataSet, .className)

    # Argument 'flavor':
    flavor <- Arguments$getCharacter(flavor, length=c(0L,1L))
  } # if (!is.null(dataSet))

  # Arguments '...':
  args <- list(...)


  this <- extend(Object(), c("AromaSeqTransform", uses("ParametersInterface")),
    .ds = dataSet,
    .tags = tags,
    .flavor = flavor,
    .args = args
  )

  setTags(this, tags)

  this
})


setMethodS3("as.character", "AromaSeqTransform", function(x, ...) {
  # To please R CMD check
  this <- x

  s <- sprintf("%s:", class(this)[1L])
  s <- c(s, sprintf("Acronym: %s", getAcronym(this)))
  s <- c(s, sprintf("Flavor: %s", getFlavor(this)))
  s <- c(s, sprintf("Number of items to process: %s", length(this)))
  s <- c(s, sprintf("Output path: %s", getPath(this, create=FALSE)))

  ds <- getInputDataSet(this)
  s <- c(s, "Input data set:")
  s <- c(s, as.character(ds))

  s <- c(s, sprintf("Parameters: %s", getParametersAsString(this)))

  GenericSummary(s)
}, protected=TRUE)



setMethodS3("length", "AromaSeqTransform", function(x) {
  length(getInputDataSet(x))
}, appendVarArgs=FALSE)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# DATA SETS
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
setMethodS3("getInputDataSet", "AromaSeqTransform", function(this, ...) {
  this$.ds
}, protected=TRUE)


setMethodS3("getOutputDataSet", "AromaSeqTransform", abstract=TRUE, protected=TRUE)

setMethodS3("findFilesTodo", "AromaSeqTransform", function(this, force=FALSE, ...) {
  # Argument 'force':
  force <- Arguments$getLogical(force)

  res <- getOutputDataSet(this, onMissing="NA")
  if (force) {
    todo <- seq_along(res)
  } else {
    isFile <- unlist(sapply(res, FUN=isFile), use.names=FALSE)
    todo <- !isFile
    todo <- which(todo)
  }
  if (length(todo) > 0L) {
    ds <- getInputDataSet(this)
    names(todo) <- getNames(ds[todo])
  }
  todo
}, protected=TRUE)


setMethodS3("isDone", "AromaSeqTransform", function(this, ...) {
  files <- findFilesTodo(this, ...)
  (length(files) == 0L)
}, protected=TRUE)


setMethodS3("getOrganism", "AromaSeqTransform", function(this, ...) {
  ds <- getInputDataSet(this)
  getOrganism(ds)
}, protected=TRUE)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# TRANSFORM SPECIFIC
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
setMethodS3("getFlavor", "AromaSeqTransform", function(this, ...) {
  this$.flavor
}, protected=TRUE)


setMethodS3("getOptionalArguments", "AromaSeqTransform", function(this, ...) {
  this$.args
}, protected=TRUE)

# TODO: Move this to aroma.core::ParametersInterface?! /HB 2013-11-16
setMethodS3("getParameters", "AromaSeqTransform", function(this, ...) {
  params <- NextMethod("getParameters")
  params <- c(params, getOptionalArguments(this, ...))
  params
}, protected=TRUE)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# NAMES, TAGS, FULL NAMES, ACRONYMS ETC.
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
setMethodS3("getAcronym", "AromaSeqTransform", function(this, case=c("upper", "lower"), ...) {
  # Argument 'case':
  case <- match.arg(case)

  # Create a default acronyme for any class from its class name
  # by extracting all capital letters and pasting them together,
  # e.g. AromaSeqTransform => AST.
  name <- class(this)[1L]
  name <- capitalize(name)
  name <- strsplit(name, split="")[[1L]]
  name <- name[(toupper(name) == name)]
  name <- paste(name, collapse="")

  if (case == "lower") {
    name <- tolower(name)
  }
}, protected=TRUE)


setMethodS3("getName", "AromaSeqTransform", function(this, ...) {
  ds <- getInputDataSet(this)
  getName(ds)
}, protected=TRUE)


setMethodS3("getAsteriskTags", "AromaSeqTransform", function(this, flavor=FALSE, collapse=NULL, ...) {
  tags <- getAcronym(this)
  if (flavor) tags <- c(tags, getFlavor(this))
  tags <- paste(tags, collapse=collapse)
  tags
}, protected=TRUE)


setMethodS3("getTags", "AromaSeqTransform", function(this, collapse=NULL, ...) {
  # "Pass down" tags from input data set
  ds <- getInputDataSet(this)
  tags <- getTags(ds, collapse=collapse)

  # Get class-specific tags
  tags <- c(tags, this$.tags)

  # Update default (aka "asterisk") tags
  asteriskTags <- getAsteriskTags(this)
  if (!is.null(asteriskTags)) {
    tags[tags == "*"] <- paste(asteriskTags, collapse=",")
  }
  tags <- tags[nchar(tags) > 0L]

  # Collapsed or split?
  if (!is.null(collapse)) {
    tags <- paste(tags, collapse=collapse)
  } else {
    tags <- unlist(strsplit(tags, split=","))
  }

  if (length(tags) == 0L) {
    tags <- NULL
  }

  tags
}, protected=TRUE)


setMethodS3("setTags", "AromaSeqTransform", function(this, tags="*", ...) {
  # Argument 'tags':
  if (!is.null(tags)) {
    tags <- Arguments$getCharacters(tags)
    tags <- trim(unlist(strsplit(tags, split=",")))
    tags <- tags[nchar(tags) > 0L]
  }
  this$.tags <- tags
  invisible(this)
}, protected=TRUE)


setMethodS3("getFullName", "AromaSeqTransform", function(this, ...) {
  name <- getName(this)
  tags <- getTags(this)
  fullname <- paste(c(name, tags), collapse=",")
  fullname <- gsub("[,]$", "", fullname)
  fullname
}, protected=TRUE)



# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# FILE PATHS
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
setMethodS3("getRootPath", "AromaSeqTransform", function(this, ...) {
  name <- getAcronym(this, case="lower")
  sprintf("%sData", name)
}, protected=TRUE)


setMethodS3("getPath", "AromaSeqTransform", function(this, create=TRUE, ...) {
  # Create the (sub-)directory tree for the data set

  # Root path
  rootPath <- getRootPath(this)

  # Full name
  fullname <- getFullName(this)

  # Organism
  organism <- getOrganism(this)

  # The full path
  path <- filePath(rootPath, fullname, organism)

  if (create) {
    path <- Arguments$getWritablePath(path)
  } else {
    path <- Arguments$getReadablePath(path, mustExist=FALSE)
  }

  # Verify that it is not the same as the input path
  ds <- getInputDataSet(this)
  inPath <- getPath(ds)
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath)
  }

  path
}, protected=TRUE)


# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
# ABSTRACT METHODS
# = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = =
setMethodS3("process", "AromaSeqTransform", abstract=TRUE)
