setMethodS3("doQDNASeq", "BamDataFile", function(df, binWidth, log=TRUE, mappability=50, blacklist=0, residual=2, bases=0, ..., force=FALSE, verbose=FALSE) {
  require("Biobase") || throw("Package not loaded: Biobase"); # combine()
  pkgName <- "qdnaseq";
  require(pkgName, character.only=TRUE) || throw("Package not loaded: qdnaseq");
  getBinAnnotations <- binReadCounts <- correctBins <- normalizeBins <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'binWidth':
  binWidth <- Arguments$getInteger(binWidth, range=c(0.1, 10e3));

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "QDNASeq");
  verbose && print(verbose, df);

  verbose && enter(verbose, "QDNASeq/Retrieve QDNASeq bin annotation");
  bins <- getBinAnnotations(binWidth);
  verbose && print(verbose, bins);
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNASeq/Reading and binning data");
  pathname <- getPathname(df);
  data <- binReadCounts(bins, bamfiles=pathname, cache=TRUE, force=force);
  verbose && print(verbose, data);
  bins <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNASeq/Correcting bin counts for GC content and mappability");
  dataC <- correctBins(data);
  verbose && print(verbose, dataC);
  data <- NULL;
  verbose && exit(verbose);

  verbose && enter(verbose, "QDNASeq/Normalization bin copy numbers");
  dataN <- normalizeBins(dataC, logTransform=log);
  verbose && print(verbose, dataN);
  dataC <- NULL;
  verbose && exit(verbose);

  verbose && exit(verbose);

  dataN;
}) # doQDNASeq()



setMethodS3("doQDNASeq", "BamDataSet", function(dataSet, binWidth, ..., verbose=FALSE) {
  pkgName <- "qdnaseq";
  require(pkgName, character.only=TRUE) || throw("Package not loaded: qdnaseq");
  getBinAnnotations <- NULL;

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'binWidth':
  binWidth <- Arguments$getInteger(binWidth);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "QDNASeq");
  verbose && print(verbose, dataSet);

  verbose && enter(verbose, "QDNASeq/Retrieve QDNASeq bin annotation");
  bins <- getBinAnnotations(binWidth);
  verbose && print(verbose, bins);
  verbose && exit(verbose);

  res <- NULL;
  for (ii in seq_along(dataSet)) {
    df <- getFile(dataSet, ii);
    verbose && enter(verbose, sprintf("Sample %d ('%s') of %d", ii, getName(df), length(dataSet)));

    dataN <- doQDNASeq(df, binWidth=binWidth, ..., verbose=less(verbose,1));
    verbose && print(verbose, dataN);

    if (is.null(res)) {
      res <- dataN;
    } else {
      res <- combine(res, dataN);
    }
    dataN <- NULL;

    verbose && exit(verbose);
  } # for (ii ...)

  verbose && print(verbose, res);

  verbose && exit(verbose);

  res;
}) # doQDNASeq()


setMethodS3("doQDNASeq", "FastqDataSet", function(ds, reference, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'reference':
  if (inherits(reference, "FastaReferenceFile")) {
  } else if (inherits(reference, "BwaIndexSet")) {
  } else {
    throw("Argument 'reference' should either be of class 'FastaReferenceFile' or 'BwaIndexSet': ", class(reference)[1L]);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "QDNASeq");

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # BWA 'aln' with options '-n 2' and '-q 40'.
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "QDNASeq/BWA alignment");

  # Retrieve BWA index set?
  if (inherits(reference, "BWAIndexSet")) {
    is <- reference;
  } else if (inherits(reference, "FastaReferenceFile")) {
    verbose && enter(verbose, "QDNASeq/BWA alignment/Retrieving index set");
    fa <- reference;
    verbose && print(verbose, fa);
    is <- buildBwaIndexSet(fa, method="is", verbose=verbose);
    verbose && print(verbose, is);
    verbose && exit(verbose);
    # Not needed anymore
    fa <- NULL;
  }
    # Not needed anymore
  reference <- NULL;

  alg <- BwaAlignment(ds, indexSet=is, n=2, q=40);
  verbose && print(verbose, alg);

  bs <- process(alg, verbose=verbose);
  verbose && print(verbose, bs);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Remove duplicated reads using Picard
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "QDNASeq/Remove duplicated reads");
  dr <- PicardDuplicateRemoval(bs);
  verbose && print(verbose, dr);

  bsU <- process(dr, verbose=verbose);
  verbose && print(verbose, bsU);

  verbose && exit(verbose);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # QDNASeq copy number estimation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  verbose && enter(verbose, "QDNASeq/copy number estimation");
  cns <- doQDNASeq(bsU, ..., verbose=verbose);
  verbose && print(verbose, cns);
  verbose && exit(verbose);

  verbose && exit(verbose);

  cns;
}) # doQDNASeq()


############################################################################
# HISTORY:
# 2013-07-11
# o Added doQDNASeq() for FastqDataSet, which leverages ditto for
#   BamDataSet.
# 2013-07-03
# o Added to verbose statements.
# 2013-07-02
# o Created.
############################################################################
