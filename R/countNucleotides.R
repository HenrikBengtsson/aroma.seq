setMethodS3("countNucleotides", "BamDataFile", function(bam, loci, ..., cache=FALSE, force=!cache, verbose=FALSE) {
  use("Rsamtools");
  use("Biostrings");

  # Need to use S4 generic as.matrix() instead of S3 one in 'base'
  ns <- getNamespace("Biostrings");
  as.matrix <- get("as.matrix", envir=ns, mode="function");


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'loci':
  loci <- Arguments$getInstanceOf(loci, "data.frame");

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Counting nucleotides in BAM file");

  pathname <- getPathname(bam);
  verbose && cat(verbose, "BAM pathname: ", pathname);
  if (!isFile(pathname)) {
    throw("BAM file does not exist: ", pathname);
  }

  verbose && cat(verbose, "Loci to inspect:");
  verbose && str(verbose, loci);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Check for memoized results
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (cache) {
    bamZ <- getChecksumFile(bam)
    key <- list(method="countNucleotides", class=class(bam)[1L], bam=readChecksum(bamZ), loci=loci)
    dirs <- c("aroma.seq", "countNucleotides", getOrganism(bam))
    counts <- loadCache(key=key, dirs=dirs);
    if (!force && !is.null(counts)) {
      verbose && cat(verbose, "Found cached results. Skipping.");
      verbose && exit(verbose);
      return(counts);
    }
  }


  # Get (chromosome, position)
  chr <- loci[,1L,drop=TRUE];
  pos <- loci[,2L,drop=TRUE];
  names <- rownames(loci);
  chr <- as.character(chr);
  if (!is.numeric(pos)) {
    throw("Second column of argument 'loci' is not numeric: ", mode(pos));
  }

  # Allocate allele counts for A, C, G, T and unknowns ("N")
  bases <- c("A", "C", "G", "T", "N")
  counts <- matrix(NA_integer_, nrow=nrow(loci), ncol=length(bases));
  colnames(counts) <- bases;
  rownames(counts) <- names;


  # SANITY CHECK: Currently only a single chromosome is supported
  uchr <- unique(chr);
  if (length(uchr) > 1L) {
    throw("Detected multiple chromosomes in argument 'loci'. Currently only a single chromosome can be scanned at the same time: ", hpaste(uchr));
  }

  chrCC <- chr[1L];
  verbose && cat(verbose, "Chromosome: ", chrCC);

  knownChromosomes <- getTargetNames(bam);
  if (!is.element(chrCC, knownChromosomes)) {
    throw(sprintf("Target sequence/chromosome '%s' is not among the known ones: %s", chrCC, hpaste(knownChromosomes)));
  }

  # Look a the loci for the chromosome of interest
  idxsCC <- which(chr == chrCC);
  posCC <- pos[idxsCC];
  namesCC <- names[idxsCC];

  # Setup RangesList for SNPs on chromosome of interest
  which <- RangesList(IRanges(start=posCC, width=1L, names=namesCC));
  names(which)[1L] <- chrCC;

  # Scan BAM file
  verbose && enter(verbose, "Calling scanBam()");
  params <- ScanBamParam(which=which, what=scanBamWhat());
  verbose && cat(verbose, "Parameters to scanBam():");
  verbose && print(verbose, params);
  res <- scanBam(pathname, param=params);
  verbose && cat(verbose, "Number of positions read: ", length(res));
  # Sanity check
  stopifnot(length(res) == length(posCC));
  verbose && exit(verbose);

  # Count alleles
  verbose && writeRaw(verbose, "Counting alleles: [0%]");
  for (jj in seq_along(res)) {
    if (verbose && (jj %% 100 == 0)) writeRaw(verbose, ".");
    resT <- res[[jj]];

    offset <- posCC[jj] - resT$pos + 1L;
    if (length(offset) > 0L) {
      seq <- resT$seq;
      seq <- as.matrix(seq);
      alleles <- rowCollapse(seq, idxs=offset);
      alleles <- factor(alleles, levels=bases);
      countsJJ <- table(alleles, dnn=NULL);
      idxJJ <- idxsCC[jj];
      counts[idxJJ,] <- countsJJ;
    }

    # Not needed anymore
    res[[jj]] <- NA;
    resT <- NULL;
  } # for (jj ...)
  res <- NULL; # Not needed anymore
  verbose && writeRaw(verbose, "[100%]\n");

  verbose && cat(verbose, "Allele counts:");
  verbose && str(verbose, counts);

  if (cache) {
    saveCache(counts, key=key, dirs=dirs);
  }

  verbose && exit(verbose);

  counts;
}) # countNucleotides()


setMethodS3("countNucleotides", "BamDataSet", function(bams, loci, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Arguments 'loci':
  loci <- Arguments$getInstanceOf(loci, "data.frame")

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Counting nucleotides in BAM set")
  verbose && print(verbose, bams)

  verbose && cat(verbose, "Loci to inspect:")
  verbose && str(verbose, loci)

  counts <- NULL

  for (ii in seq_along(bams)) {
    bam <- bams[[ii]]
    verbose && enterf(verbose, "File #%d ('%s') of %d", ii, getName(bam), length(bams))
    verbose && print(verbose, bam)

    countsII <- countNucleotides(bam, loci=loci, ..., verbose=less(verbose, 5))
    dimII <- dim(countsII)

    if (is.null(counts)) {
      counts <- matrix(0L, nrow=dimII[1L], ncol=dimII[2L])
      dimnames(counts) <- dimnames(countsII)
    }
    stopifnot(identical(dimII, dim(counts)))

    countsII[is.na(countsII)] <- 0L
    counts <- counts + countsII

    countsII <- NULL # Not needed anymore

    verbose && exit(verbose)
  } # for (ii ...)

  total <- rowSums(counts, na.rm=TRUE)
  counts[(total == 0L),] <- NA_integer_

  verbose && cat(verbose, "Allele counts:")
  verbose && str(verbose, counts)

  verbose && exit(verbose)

  counts;
}) # countNucleotides()


############################################################################
# HISTORY:
# 2014-09-28
# o Added countNucleotides() for BamDataSet.
# o Added support for memoization to countNucleotides() for BamDataFile.
# 2014-06-14
# o Added to aroma.seq.
# o Created from BAF,chr{{chr}}.md.rsp.
############################################################################
