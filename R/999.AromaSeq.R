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
# @author
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
#  \item{...}{Not used.}
# }
#
# \value{
#   Returns a @logical named @character @vector.
# }
#
# \examples{
#   # Display which sequencing tools are supported by the package
#   print(capabilitiesOf(aroma.seq))
#
#   # Check whether BWA is supported
#   print(isCapableOf(aroma.seq, "bwa"))
# }
#
# @author
#*/###########################################################################
setMethodS3("capabilitiesOf", "AromaSeq", function(static, what=NULL, ...) {
  res <- list();

  # General software frameworks
  res$java <- !is.null(findJava(mustExists=FALSE));
  res$perl <- !is.null(findPerl(mustExists=FALSE));

  # Sequencing tools
  res$bowtie2 <- !is.null(findBowtie2(mustExists=FALSE));
  res$bwa <- !is.null(findBWA(mustExists=FALSE));
  res$gatk <- !is.null(findGATK(mustExists=FALSE));
  res$picard <- !is.null(findPicard(mustExists=FALSE));
  res$samtools <- !is.null(findSamtools(mustExists=FALSE));
  res$tophat1 <- !is.null(findTopHat(version=1, mustExists=FALSE));
  res$tophat2 <- !is.null(findTopHat(version=2, mustExists=FALSE));

  # Order lexicographically
  o <- order(names(res));
  res <- res[o];

  # Coerce into a named character vector
  res <- unlist(res);

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



############################################################################
# HISTORY:
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
