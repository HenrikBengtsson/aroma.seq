###########################################################################/**
# @RdocDefault bwaSampe
#
# @title "Generates BWA-backtrack paired-end (PE) alignments via 'bwa sampe'"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathnameSAI}{A @character @vector of two SAI files.}
#   \item{pathnameFQ}{A @character @vector of two FASTQ or BAM files.}
#   \item{indexPrefix}{The pathname prefix to the BWA index files.}
#   \item{pathnameD}{The destination pathname.}
#   \item{...}{Additional arguments specifying BWA 'sampe' switches
#     passed to @see "systemBWA".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \examples{\dontrun{
#   pathnameFA <- "annotationData/organisms/Lambda_phage/lambda_virus.fa"
#   bwaIndex(pathnameFA)
#
#   pathnameSAI <- "bwaData/LambdaVirusExample/Lambda_phage/reads_1.sai"
#   pathnameFQ <- "fastqData/LambdaVirusExample/Lambda_phage/reads_1.fq"
#   pathnameD <- "bwaData/LambdaVirusExample/Lambda_phage/reads_1.sam"
#   bwaSampe(pathnameSAI=pathnameSAI, pathnameFQ=pathnameFQ,
#            pathnameFA=pathnameFA, pathnameD=pathnameD)
# }}
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setMethodS3("bwaSampe", "default", function(pathnameSAI, pathnameFQ, indexPrefix, pathnameD, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathnameSAI':
  pathnameSAI <- Arguments$getReadablePathnames(pathnameSAI, length=c(2,2))
  .stop_if_not(pathnameSAI[2L] != pathnameSAI[1L])

  # Argument 'pathnameFQ':
  pathnameFQ <- Arguments$getReadablePathnames(pathnameFQ, length=c(2,2))
  .stop_if_not(pathnameFQ[2L] != pathnameFQ[1L])

  # Argument 'indexPrefix':
  dummy <- Arguments$getReadablePath(getParent(indexPrefix))

  # Argument 'pathnameD':
  pathnameD <- Arguments$getWritablePathname(pathnameD)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Running BWA 'sampe'")

  # Assert that input files are not overwritten
  .stop_if_not(!getAbsolutePath(pathnameD) %in% getAbsolutePath(pathnameSAI))
  .stop_if_not(!getAbsolutePath(pathnameD) %in% getAbsolutePath(pathnameFQ))
##  .stop_if_not(getAbsolutePath(pathnameD) != getAbsolutePath(pathnameFA))

  res <- systemBWA("sampe", "f"=shQuote(pathnameD), shQuote(indexPrefix), shQuote(pathnameSAI), shQuote(pathnameFQ), ..., verbose=less(verbose, 10))

  verbose && exit(verbose)

  res
}) # bwaSampe()
