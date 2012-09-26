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
  res$bowtie2 <- !is.null(findBowtie2(mustExists=FALSE));
  res$bwa <- !is.null(findBWA(mustExists=FALSE));
  res$samtools <- !is.null(findSamtools(mustExists=FALSE));
  res <- unlist(res);

  if (!is.null(what)) {
    res <- res[what];
  }

  res;
}, static=TRUE)


setMethodS3("isCapableOf", "AromaSeq", function(static, what, ...) {
  capabilitiesOf(static, what=what, ...);
})



############################################################################
# HISTORY:
# 2012-09-25
# o Added capabilitiesOf().
# o Added 'samtools' to AromaSeq$isCapableOf().
# 2012-09-24
# o Created.
############################################################################
