###########################################################################/**
# @RdocClass HTSeqCountDataFile
#
# @title "The HTSeqCountDataFile class"
#
# \description{
#  @classhierarchy
#
#  A HTSeqCountDataFile represents a HT-Seq [1] (htseq-count) output data file.
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::TabularTextFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# @author "HB")
#
# \references{
#  [1] Simon Anders, \emph{HTSeq: Analysing high-throughput sequencing
#      data with Python}, EMBL, Jan 2014.
#      \url{http://www-huber.embl.de/users/anders/HTSeq/} \cr
# }
#*/###########################################################################
setConstructorS3("HTSeqCountDataFile", function(...) {
  extend(TabularTextFile(..., columnNames=c("transcript", "count"), .verify=FALSE), c("HTSeqCountDataFile", uses("AromaSeqDataFile")))
})


##############################################################################
# HISTORY:
# 2014-01-24
# o Created.
##############################################################################
