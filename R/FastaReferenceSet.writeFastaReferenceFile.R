###########################################################################/**
# @set class=FastaReferenceSet
# @RdocMethod writeFastaReferenceFile
# @alias writeFastaReferenceFile
#
# @title "Writes the content of multiple FASTA files into a new FASTA file"
#
# \description{
#   @get "title".  Existing files are left unmodified.
# }
#
# @synopsis
#
# \arguments{
#  \item{filename, path}{The filename and path of the generated
#        @see "FastaReferenceFile".}
#  \item{...}{Not used.}
#  \item{overwrite}{If @TRUE and the output file already exists, then
#    it is overwritten, otherwise an exception is thrown.}
#  \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# \value{
#   Returns a @see "FastaReferenceFile".
# }
#
# @author "HB"
#
# @keyword IO
# @keyword internal
#*/###########################################################################
setMethodS3("writeFastaReferenceFile", "FastaReferenceSet", function(this, filename=sprintf("%s,%s.fa", getOrganism(this), getChecksum(this)), path=file.path("annotationData", "organisms", getOrganism(this)), ..., overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'filename' and 'path':
  pathname <- Arguments$getWritablePathname(filename, path=path,
                                            mustNotExist=!overwrite)

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enterf(verbose, "Writing %s to one FastaReferenceFile", class(this)[1L])
  verbose && cat(verbose, "Input data set:")
  verbose && print(verbose, this)
  verbose && cat(verbose, "Output pathname: ", pathname)

  # Never overwrite one of the input FASTA files
  if (overwrite && isFile(pathname)) {
    fa <- FastaReferenceFile(pathname)
    for (ii in seq_along(this)) {
      if (equals(this[[ii]], fa)) {
        throw("Cannot write FASTA file. The output pathname refers to one of the input files: ", pathname)
      }
    }
  }

  # Size of buffer when appending a file in chunks
  BFRSIZE <- 100e6

  # Write to a temporary file
  pathnameT <- pushTemporaryFile(pathname, verbose=verbose)

  # Output connection
  con <- file(pathnameT, open="w+b")
  on.exit({
    if (!is.null(con)) close(con)
  }, add=TRUE)

  # Total number of bytes written
  total <- 0

  for (ii in seq_along(this)) {
    fa <- this[[ii]]
    verbose && enterf(verbose, "File #%d ('%s') of %d", ii, getName(fa), length(this))
    verbose && print(verbose, fa)

    # Append current file in chunks
    pathnameFA <- getPathname(fa)
    conII <- gzfile(pathnameFA, open="rb")
    appendNewline <- FALSE
    totalII <- 0L
    n <- Inf
    repeat {
      bfr <- readBin(conII, what="raw", n=BFRSIZE)
      n <- length(bfr)
      if (n == 0L) break
      writeBin(bfr, con=con)
      totalII <- totalII + n
      appendNewline <- (bfr[n] != charToRaw("\n"))
      bfr <- NULL
    }
    close(conII); conII <- NULL

    # Total number of bytes written
    total <- total + totalII
    verbose && printf(verbose, "Total number of bytes: %.0f (added %.0f)\n", total, totalII)

    # Append missing newline?
    if (appendNewline) {
      writeBin("\n", con=con, size=1L)
      total <- total + 1
    }

    verbose && exit(verbose)
  } # for (ii ...)

  # Close output connection
  close(con); con <- NULL

  # Renaming temporary file
  pathname <- popTemporaryFile(pathnameT, verbose=verbose)

  # Return result
  res <- FastaReferenceFile(pathname)
  verbose && print(verbose, res)

  # Make sure to checksum file exists
  resZ <- getChecksumFile(res)
  verbose && print(verbose, resZ)

  verbose && exit(verbose)

  res
}) # writeFastaReferenceFile()
