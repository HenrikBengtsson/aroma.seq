###########################################################################/**
# @RdocClass SraDataFile
#
# @title "The abstract SraDataFile class"
#
# \description{
#  @classhierarchy
#
#  A SraDataFile object represents a Sequence Read Archive (SRA) file [1].
# }
#
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "R.filesets::GenericDataFile".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \references{
#  [1] \emph{The Sequence Read Archive (SRA)},
#      \url{http://www.ncbi.nlm.nih.gov/sra/}, 2014.\cr
# }
#
# \seealso{
#   An object of this class is typically part of an
#   @see "SraDataSet".
# }
#
# @author "HB"
#*/###########################################################################
setConstructorS3("SraDataFile", function(...) {
  extend(GenericDataFile(...), c("SraDataFile", uses("AromaSeqDataFile")))
})


setMethodS3("as.character", "SraDataFile", function(x, ...) {
  this <- x
  s <- NextMethod("as.character")
  s
}, protected=TRUE)


setMethodS3("fastqDump", "SraDataFile", function(this, path=".", ..., skip=TRUE, overwrite=!skip, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'path':
  path <- Arguments$getWritablePath(path)

  # Argument 'skip':
  skip <- Arguments$getLogical(skip)

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Writing SRA file as FASTQ file")
  verbose && cat(verbose, "Input file:")
  verbose && print(verbose, this)
  pathnameSRA <- getPathname(this)
  pathnameSRA <- Arguments$getReadablePathname(pathnameSRA)

  fullname <- getFullName(this)
  filenameFQ <- sprintf("%s.fastq", fullname)
  pathnameFQ <- Arguments$getWritablePathname(filenameFQ, path=path, mustNotExist=FALSE)
  verbose && cat(verbose, "Output pathname: ", pathnameFQ)

  pathnameSRA <- getAbsolutePath(pathnameSRA)
  pathnameFQ <- getAbsolutePath(pathnameFQ)

  if (isFile(pathnameFQ)) {
    if (skip) {
      verbose && cat(verbose, "Detected existing FASTQ file. Skipping.")
    } else {
      if (overwrite) {
        verbose && cat(verbose, "Deleting existing FASTQ file (overwrite=TRUE): ", pathnameFQ)
        file.remove(pathnameFQ)
      } else {
        throw("Output FASTQ file already exists: ", pathnameFQ)
      }
    }
  } else {
    skip <- FALSE
  }

  if (!skip) {
    verbose && enter(verbose, "Calling fastq-dump of the SRA Toolkit")
    
    pathT <- file.path(path, sprintf("%s.fastq.tmp", fullname))
    pathT <- Arguments$getWritablePath(pathT)

    opwd <- setwd(pathT)
    on.exit({
      setwd(opwd)
    }, add=TRUE)

    res <- systemSraToolkit("fastq-dump", pathnameSRA, ..., verbose=verbose)

    verbose && cat(verbose, "Result code:")
    verbose && str(verbose, res)
    verbose && exit(verbose)

    # Assert output
    pathnameFQT <- Arguments$getReadablePathname(filenameFQ)

    # Renam temporary path
    renameFile(pathnameFQT, pathnameFQ, overwrite=TRUE)
    
    setwd(opwd)

    removeDirectory(pathT)
  }

  fq <- FastqDataFile(pathnameFQ)
  verbose && print(verbose, fq)

  verbose && exit(verbose)

  fq
})
