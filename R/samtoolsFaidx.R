###########################################################################/**
# @RdocDefault samtoolsFaidx
#
# @title "Calls the samtools 'faidx' command"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{pathname}{A reference fasta file.}
#   \item{...}{Additional arguments specifying samtools 'faidx' switches
#     passed to @see "systemSamtools".}
#   \item{verbose}{See @see "R.utils::Verbose".}
# }
#
# @author "TT"
#
# @keyword internal
#*/###########################################################################
setMethodS3("samtoolsFaidx", "default", function(pathname, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'pathname':
  pathname <- Arguments$getReadablePathname(pathname);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Running samtools 'faidx'");

  res <- systemSamtools("faidx", shQuote(pathname), ..., verbose=less(verbose, 10));

  verbose && exit(verbose);

  res;
}) # samtoolsFaidx()


# From http://samtools.sourceforge.net/samtools.shtml:
# faidx 	samtools faidx <ref.fasta> [region1 [...]]
# Index reference sequence in the FASTA format or extract subsequence from indexed reference sequence. If no region is specified, faidx will index the file and create <ref.fasta>.fai on the disk. If regions are speficified, the subsequences will be retrieved and printed to stdout in the FASTA format. The input file can be compressed in the RAZF format.


############################################################################
# HISTORY:
# 2014-03-10 [HB]
# o ROBUSTNESS: Now samtoolsFaidx() uses shQuote() for all pathnames.
# 2013-11-15
# o Created from samtoolsView.R
############################################################################
