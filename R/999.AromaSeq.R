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
  res <- static$.capabilities;
  if (force || is.null(res)) {
    res <- list();

    # General software frameworks
    res$java <- !is.null(findJava(mustExist=FALSE));
    res$perl <- !is.null(findPerl(mustExist=FALSE));
    res$python <- !is.null(findPython(mustExist=FALSE));

    # Sequencing tools
    res$bowtie2 <- !is.null(findBowtie2(mustExist=FALSE));
    res$bwa <- !is.null(findBWA(mustExist=FALSE));
    res$gatk <- !is.null(findGATK(mustExist=FALSE));
    res$picard <- !is.null(findPicard(mustExist=FALSE));
    res$fastqc <- !is.null(findFastQC(mustExist=FALSE));
    res$samtools <- !is.null(findSamtools(mustExist=FALSE));
    res$tophat1 <- !is.null(findTopHat1(mustExist=FALSE));
    res$tophat2 <- !is.null(findTopHat2(mustExist=FALSE));
    res$htseq <- !is.null(findHTSeq(mustExist=FALSE));

    # Order lexicographically
    o <- order(names(res));
    res <- res[o];

    # Coerce into a named character vector
    res <- unlist(res);

    # Record
    static$.capabilities <- res;
  }

  if (!is.null(what)) {
    res <- res[what];
  }

  res;
}, static=TRUE)


setMethodS3("isCapableOf", "AromaSeq", function(static, what, ...) {
  capabilitiesOf(static, what=what, ...);
})


setMethodS3("setupTests", "AromaSeq", function(static, path="redundancyTests/", ...) {
  # Argument 'path':
  path <- Arguments$getWritablePath(path);

  # Get the setup script
  pathT <- system.file("testScripts", "setup", package=getName(static));
  pathname <- Arguments$getReadablePathname("00a.setup.R", path=pathT);

  opwd <- getwd();
  setwd(path);
  on.exit(setwd(opwd));

  # Setup test directory
  source(pathname);

  path;
})

# \references{
#   \url{http://bioconductor.org/packages/release/BiocViews.html#___AnnotationData}
# }
setMethodS3("getKnownOrganisms", "AromaSeq", function(static, ...) {
  c(
    "DrosophilaMelanogaster",
    "EscherichiaColi",
    "HomoSapiens",
    "LambdaPhage",
    "MusMusculus"
  );
}, protected=TRUE)


setMethodS3("getOrganism", "Arguments", function(static, organism, ...) {
  # Argument 'organism':
  organism <- Arguments$getCharacter(organism, length=c(1L,1L));
  knownOrganisms <- getKnownOrganisms(aroma.seq);
  unknown <- organism[!is.element(organism, knownOrganisms)];
  if (length(unknown) > 0L) {
    throw("Unknown organism: ", organism);
  }

  organism;
}, protected=TRUE)


setMethodS3("skeleton", "AromaSeq", function(static, dataSet="MyDatSet", organism="HomoSapiens", ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'dataSet':
  dataSet <- Arguments$getCharacter(dataSet);

  # Argument 'organism':
  knownOrganisms <- getKnownOrganisms(static);
  organism <- Arguments$getOrganism(organism);

  if (dataSet == organism) {
    warning("Did you really mean to name the data set the same as the organism?: ", dataSet);
  }

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # annotationData/organisms/<organism>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- file.path("annotationData", "organisms", organism);
  path <- Arguments$getWritablePath(path);
  pathname <- file.path(path, "README.txt");
  if (!isFile(pathname)) {
    cat("Copy or link to the FASTA reference file in this directory.\n", file=pathname);
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # fastqData/<DataSet>/<organism>/
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  path <- file.path("fastqData", dataSet, organism);
  path <- Arguments$getWritablePath(path);
  pathname <- file.path(path, "README.txt");
  if (!isFile(pathname)) {
    cat("Copy or link to the FASTQ read files in this directory.\n", file=pathname);
  }

  invisible(TRUE);
}) # skeleton()


############################################################################
# HISTORY:
# 2014-02-28
# o Added 'fastqc' to capabilities.
# 2013-11-08
# o Added skeleton() for AromaSeq, e.g.
#   skeleon(aroma.seq, "MyDataSet", "HomoSapiens").
# 2013-10-30
# o Added 'python' to capabilities.
# 2013-07-19
# o SPEEDUP: Now the results of capabilitiesOf() are cached.
# 2013-07-03
# o Added 'HTSeq' to capabilitiesOf().
# 2013-04-01
# o Added 'tophat1' and 'tophat2' to capabilitiesOf().
# 2012-09-27
# o Added setupTests() for AromaSeq.
# 2012-09-25
# o Added capabilitiesOf().
# o Added 'samtools' to AromaSeq$isCapableOf().
# 2012-09-24
# o Created.
############################################################################
