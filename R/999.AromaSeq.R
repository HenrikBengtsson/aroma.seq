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
# @RdocMethod isCapableOf
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
#   # Display which sequencing tools are supported
#   print(AromaSeq$isCapableOf())
#
#   # Check whether BWA is supported
#   print(AromaSeq$isCapableOf("bwa"))
# }
#
# @author
#*/###########################################################################
setMethodS3("isCapableOf", "AromaSeq", function(static, what=NULL, ...) {
  res <- list();
  res$bowtie2 <- !is.null(findBowtie2(mustExists=FALSE));
  res$bwa <- !is.null(findBWA(mustExists=FALSE));
  res <- unlist(res);

  if (!is.null(what)) {
    res <- res[what];
  }

  res;
}, static=TRUE)



############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################
