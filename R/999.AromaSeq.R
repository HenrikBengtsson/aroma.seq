###########################################################################/**
# @RdocClass AromaSeq
#
# @title "The AromaSeq Package class"
#
# \description{
#  @classhierarchy
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Not used.}
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
setConstructorS3("AromaSeq", function(...) {
  extend(AromaPackage("aroma.seq", ...), "AromaSeq");
})



###########################################################################/**
# @RdocMethod capabilitiesOf
# @aliasmethod isCapableOf
#
# @title "Checks which tools are supported"
#
# \description{
#   @get "title".
# }
#
# @synopsis
#
# \arguments{
#  \item{what}{Optional @character @vector of which tools to check.}
#  \item{force}{If @TRUE, cached results are ignored, otherwise not.}
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @logical named @character @vector.
# }
#
# \examples{
#   # Display which tools are supported by the package
#   print(capabilitiesOf(aroma.seq))
#
#   # Check whether BWA is supported
#   print(isCapableOf(aroma.seq, "bwa"))
# }
#
# @author "HB"
#
#*/###########################################################################
setMethodS3("capabilitiesOf", "AromaSeq", function(static, what=NULL, force=FALSE, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  supports <- function(fcn, ...) {
    tryCatch({
      !is.null(fcn(mustExist=FALSE))
    }, error = function(ex) FALSE)
  } # supports()

  res <- static$.capabilities
  if (force || is.null(res)) {
    res <- list()

    # General software frameworks
    res$java <- supports(findJava)
    res$perl <- supports(findPerl)
    res$python <- supports(findPython)

    # Sequencing tools
    res$bowtie2 <- supports(findBowtie2)
    res$bwa <- supports(findBWA)
    res$gatk <- supports(findGATK)
    res$CNVkit <- supports(findCNVkit)
    res$picard <- supports(findPicard)
    res$fastqc <- supports(findFastQC)
    res$fastqDump <- supports(findFastqDump)
    res$samtools <- supports(findSamtools)
    res$sratoolkit <- supports(findSraToolkit)
    res$tophat1 <- supports(findTopHat1)
    res$tophat2 <- supports(findTopHat2)
    res$htseq <- supports(findHTSeq)

    # Order lexicographically
    withLocale({
      o <- order(names(res))
      res <- res[o]
    }, category="LC_COLLATE", locale="C")

    # Coerce into a named character vector
    res <- unlist(res)

    # Record
    static$.capabilities <- res
  }

  if (!is.null(what)) {
    res <- res[what]
  }

  res
}, static=TRUE)


setMethodS3("isCapableOf", "AromaSeq", function(static, what, ...) {
  capabilitiesOf(static, what=what, ...)
})


setMethodS3("setupTests", "AromaSeq", function(static, path="redundancyTests/", ...) {
  # Argument 'path':
  path <- Arguments$getWritablePath(path)

  # Get the setup script
  pathT <- system.file("testScripts", "setup", package=getName(static))
  pathname <- Arguments$getReadablePathname("00a.setup.R", path=pathT)

  opwd <- getwd()
  setwd(path)
  on.exit(setwd(opwd))

  # Setup test directory
  source(pathname)

  path
})

# \references{
#   \url{http://bioconductor.org/packages/release/BiocViews.html#___AnnotationData}
# }
setMethodS3("getKnownOrganisms", "AromaSeq", function(static, ...) {
  c(
    "Drosophila_melanogaster",
    "Escherichia_coli",
    "Homo_sapiens",
    "Lambda_phage",
    "Mus_musculus"
  )
}, protected=TRUE)


setMethodS3("getOrganism", "Arguments", function(static, organism, mustBeKnown=FALSE, ...) {
  # If an (aroma.*) file or a data set it passed, get the organism
  # string for that object.
  if (inherits(organism, "GenericDataFile") ||
      inherits(organism, "GenericDataFileSet")) {
    organism <- getOrganism(organism)
  }

  # Assert a single character string
  organism <- Arguments$getCharacter(organism, length=c(1L,1L))

  # Assert a known organism?
  if (mustBeKnown) {
    knownOrganisms <- getKnownOrganisms(aroma.seq)
    unknown <- organism[!is.element(organism, knownOrganisms)]
    if (length(unknown) > 0L) {
      throw("Unknown organism: ", organism)
    }
  }

  organism
}, protected=TRUE)


setMethodS3("skeleton", "AromaSeq", function(static, dataSet="MyDatSet", organism="Homo_sapiens", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet)

  # Argument 'organism':
  organism <- Arguments$getOrganism(organism, ...)

  if (dataSet == organism) {
    warning("Did you really mean to name the data set the same as the organism?: ", dataSet)
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # annotationData/organisms/<organism>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- file.path("annotationData", "organisms", organism)
  path <- Arguments$getWritablePath(path)
  pathname <- file.path(path, "README.txt")
  if (!isFile(pathname)) {
    cat("Copy or link to the FASTA reference file in this directory.\n", file=pathname)
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # fastqData/<DataSet>/<organism>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- file.path("fastqData", dataSet, organism)
  path <- Arguments$getWritablePath(path)
  pathname <- file.path(path, "README.txt")
  if (!isFile(pathname)) {
    cat("Copy or link to the FASTQ read files in this directory.\n", file=pathname)
  }

  invisible(TRUE)
}) # skeleton()
