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
#  \item{dataSet}{An @see "AromaSeqDataFileSet".}
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
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setConstructorS3("AbstractAlignment", function(dataSet=NULL, indexSet=NULL, rgSet=NULL, ...) {
  # Validate arguments
  if (!is.null(dataSet)) {
    # Argument 'dataSet':
    dataSet <- Arguments$getInstanceOf(dataSet, "AromaSeqDataFileSet");

    # Argument 'indexSet':
    indexSet <- Arguments$getInstanceOf(indexSet, "AbstractIndexSet");

    # Argument 'rgSet':
    if (is.null(rgSet)) {
      if (inherits(dataSet, "FastqDataSet")) {
        rgSet <- getSamReadGroup(dataSet)
      }
    } else {
      rgSet <- Arguments$getInstanceOf(rgSet, "SamReadGroup");
    }
  } # if (!is.null(dataSet))


  extend(AromaSeqTransform(dataSet=dataSet, ...), "AbstractAlignment",
    .indexSet = indexSet,
    .rgSet = rgSet
  )
})


setMethodS3("as.character", "AbstractAlignment", function(x, ...) {
  s <- NextMethod("as.character");
  s <- c(s, sprintf("Paired alignment: %s", isPaired(x)));
  s;
}, protected=TRUE)

setMethodS3("isPaired", "AbstractAlignment", function(this, ...) {
  ds <- getInputDataSet(this);
  isPaired(ds, ...);
}, protected=TRUE)


setMethodS3("getIndexSet", "AbstractAlignment", function(this, ...) {
  this$.indexSet;
}, protected=TRUE)

# This methods use lower case by default, e.g. bwa and bowtie2.
setMethodS3("getAcronym", "AbstractAlignment", function(this, case=c("lower", "upper"), ...) {
  # Argument 'case':
  case <- match.arg(case);

  name <- class(this)[1L];
  name <- gsub("Alignment", "", name, fixed=TRUE);
  name <- toupper(name);

  if (case == "lower") {
    name <- tolower(name);
  }
  name;
}, protected=TRUE)


setMethodS3("getAsteriskTags", "AbstractAlignment", function(this, collapse=NULL, ...) {
  tags <- NextMethod("getAsteriskTags");

  # Tags when paired-end reads are used
  if (isPaired(this)) tags <- c(tags, "pe");

  # Tags for the index set
  is <- getIndexSet(this);
  tags <- c(tags, getTags(is, collapse=NULL));

  tags <- unique(tags);
  paste(tags, collapse=collapse);
}, protected=TRUE)


setMethodS3("getOrganism", "AbstractAlignment", function(this, ...) {
  is <- getIndexSet(this);
  getOrganism(is);
}, protected=TRUE)

setMethodS3("getRootPath", "AbstractAlignment", function(this, ...) {
  "bamData";
}, protected=TRUE)
