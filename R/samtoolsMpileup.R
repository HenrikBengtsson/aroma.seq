###########################################################################/**
# @RdocDefault samtoolsMpileup
#
# @title "Calls the samtools 'mpileup' command"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{refFile}{input reference file}
#   \item{bamFile}{input BAM file (currently only one file is supported)}
#   \item{pathnameD}{output file (default 'mpileup.out')}
#   \item{...}{Additional arguments specifying samtools 'mpileup' switches
#     passed to @see "systemSamtools".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("samtoolsMpileup", "default", function(refFile, bamFile, pathnameD="mpileup.out", ..., verbose=FALSE) {

  # Support a call like this:  system(paste("samtools mpileup -uf", RefFile, BamFile, ">tmp1.out"))
  # NB from samtools mpileup help:  'Assuming diploid individuals.'  

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  refFile <- Arguments$getReadablePathname(refFile);
  bamFile <- Arguments$getReadablePathname(bamFile);
  
  # Argument 'pathnameD':
  pathnameD <- Arguments$getWritablePathname(pathnameD);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running samtools 'mpileup'");

  # Assert that input files are not overwritten
  stopifnot((getAbsolutePath(pathnameD) != getAbsolutePath(refFile)) &&
              (getAbsolutePath(pathnameD) != getAbsolutePath(bamFile)));
  
  res <- systemSamtools("mpileup", ..., refFile, bamFile, "stdout"=pathnameD, verbose=less(verbose, 10));

  verbose && exit(verbose);

  res;
}) # samtoolsMpileup()


############################################################################
# HISTORY:
# 2013-10-31
# o Created as copy of samtoolsView.R [TT]
############################################################################
