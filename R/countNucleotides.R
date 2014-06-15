setMethodS3("countNucleotides", "BamDataFile", function(bam, loci, ..., verbose=FALSE) {
  use("Rsamtools");

  # Arguments 'loci':
  loci <- Arguments$getInstanceOf(loci, "data.frame");

  pathnameBAM <- getPathname(bam);
  if (!isFile(pathnameBAM)) {
    throw("BAM file does not exist: ", pathnameBAM);
  }

  # Get (chromosome, position)
  chr <- loci[,1L,drop=TRUE];
  pos <- loci[,2L,drop=TRUE];
  names <- rownames(loci);

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

  # Look a the loci for the chromosome of interest
  idxsCC <- which(chr == chrCC);
  posCC <- pos[idxsCC];
  namesCC <- names[idxsCC];

  # Setup RangesList for SNPs on chromosome of interest
  which <- RangesList(IRanges(start=posCC, width=1L, names=namesCC));
  names(which)[1L] <- sprintf("%d", chrCC);

  # Scan BAM file
  params <- ScanBamParam(which=which, what=scanBamWhat());
  res <- scanBam(pathnameBAM, param=params);

  # Count alleles
  for (jj in seq_along(res)) {
    resT <- res[[jj]];
    offset <- posCC[jj] - resT$pos + 1L;
    if (length(offset) > 0L) {
      seq <- as.matrix(resT$seq);
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

  counts;
}) # countNucleotides()


############################################################################
# HISTORY:
# 2014-06-14
# o Added to aroma.seq.
# o Created from BAF,chr{{chr}}.md.rsp.
############################################################################
