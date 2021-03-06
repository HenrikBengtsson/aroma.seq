## ac <- GatkAlleleCounting(bs, targetUgp=ugp, fa=fa)
## print(ac)
## dsC <- process(ac, verbose=verbose)
## print(dsC)

setConstructorS3("GatkAlleleCounting", function(dataSet=NULL, targetUgp=NULL, fa=NULL, dropEmpty=TRUE, ..., .reqSetClass="BamDataSet") {
  use("aroma.cn")

  # Argument 'targetUgp':
  if (!is.null(targetUgp)) {
    targetUgp <- Arguments$getInstanceOf(targetUgp, "AromaUgpFile")
  }

  # Argument 'fa':
  if (!is.null(fa)) {
    fa <- Arguments$getInstanceOf(fa, "FastaReferenceFile")
  }

  # Argument 'dropEmpty':
  dropEmpty <- Arguments$getLogical(dropEmpty)


  extend(AromaTransform(dataSet=dataSet, ..., .reqSetClass=.reqSetClass), "GatkAlleleCounting",
    .targetUgp = targetUgp,
    .fa = fa,
    .dropEmpty = dropEmpty
  )
}) # GatkAlleleCounting()


setMethodS3("getParameters", "GatkAlleleCounting", function(this, ...) {
  params <- list(
    targetUgp = this$.targetUgp,
    fa = this$.fa,
    dropEmpty = this$.dropEmpty
  )
  params
}, protected=TRUE)


setMethodS3("getRootPath", "GatkAlleleCounting", function(this, ...) {
  "gatkData"
}, protected=TRUE)


setMethodS3("getPath", "GatkAlleleCounting", function(this, ...) {
  path <- NextMethod("getPath", create=FALSE)
  path <- dirname(path)
  params <- getParameters(this)
  targetUgp <- params$targetUgp
  chipType <- getChipType(targetUgp, fullname=FALSE)

  # The full path
  path <- filePath(path, chipType)
  path <- Arguments$getWritablePath(path)

  # Verify that it is not the same as the input path
  inPath <- getPath(getInputDataSet(this))
  if (getAbsolutePath(path) == getAbsolutePath(inPath)) {
    throw("The generated output data path equals the input data path: ", path, " == ", inPath)
  }

  path
}, protected=TRUE)

setMethodS3("getExpectedOutputFullnames", "GatkAlleleCounting", function(this, ...) {
  names <- NextMethod("getExpectedOutputFullNames")
  names <- paste(names, "alleleCounts", sep=",")
  names
}, protected=TRUE)


setMethodS3("getOutputDataSet0", "GatkAlleleCounting", function(this, ...) {
  NextMethod("getOutputDataSet0", className="TabularTextFileSet", pattern=".*,allel[e]*Counts[.]txt$")
}, protected=TRUE)


setMethodS3("process", "GatkAlleleCounting", function(this, ..., overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }


  verbose && enter(verbose, "Counting alleles for known SNPs")

  bams <- getInputDataSet(this)
  verbose && print(verbose, bams)

  params <- getParameters(this)
  targetUgp <- params$targetUgp
  fa <- params$fa
  dropEmpty <- params$dropEmpty

  verbose && print(verbose, targetUgp)
  verbose && print(verbose, fa)

  # Get output path
  pathD <- getPath(this)
  verbose && cat(verbose, "Output path: ", pathD)

  bedf <- NULL

  for (ii in seq_along(bams)) {
    bam <- bams[[ii]]
    verbose && enter(verbose, sprintf("Sample #%d ('%s') of %d", ii, getName(bam), length(bams)))

    filename <- sprintf("%s,alleleCounts.txt", getFullName(bam))
    pathnameD <- Arguments$getReadablePathname(filename, path=pathD, mustExist=FALSE)

    if (overwrite || !isFile(pathnameD)) {
      pathnameD <- Arguments$getWritablePathname(pathnameD, mustNotExist=!overwrite)

      # (a) Instead of having GATK build a missing FAI index file, build it here
      buildIndex(fa, verbose=verbose)

      # (b) Get the BED file representation of UGP file
      if (is.null(bedf)) {
        verbose && enter(verbose, "Writing BED file representation of UGP file")
        bedf <- writeBedDataFile(targetUgp, chrMap=c(X=23, Y=24, MT=25), verbose=verbose)
        verbose && print(verbose, bedf)
        verbose && exit(verbose)
      }

      # (c) Call GATK
      verbose && enter(verbose, "Calling GATK DepthOfCoverage")
      pathnameDT <- pushTemporaryFile(pathnameD)
      verbose && cat(verbose, "Writing to temporary file: ", pathnameDT)
      # -T / --analysis_type: Name of the tool to run
      # -I / --input_file: Pathname to BAM file
      # -R / --reference_sequence: Pathname to FASTA reference file
      # -L / --intervals: Pathname to BED file
      # --omitIntervalStatistics: Do not calculate per-interval statistics
      #   [Disabling the tabulation of interval statistics (mean, median, quartiles AND # intervals by sample by coverage) should speed up processing. This option is required in order to use -nt parallelism.]
      # --omitLocusTable: Do not calculate per-sample per-depth counts of loci
      #   [Disabling the tabulation of locus statistics (# loci covered by sample by coverage) should speed up processing.]
      # --omitPerSampleStats: Do not output the summary files per-sample
      #   [This option simply disables writing separate files for per-sample summary statistics (total, mean, median, quartile coverage per sample). These statistics are still calculated internally, so enabling this option will not improve runtime.]
      # --printBaseCounts: Add base counts to per-locus output
###      res <- systemGATK(T="DepthOfCoverage", I=getPathname(bam), R=getPathname(fa), L=getPathname(bedf), "--omitIntervalStatistics", "--omitLocusTable", "--omitPerSampleStats", "--printBaseCounts", "o"=pathnameDT, verbose=verbose)

      res <- gatk(analysisType="DepthOfCoverage", pathnameI=getPathname(bam), pathnameR=getPathname(fa), pathnameL=getPathname(bedf), "--omitIntervalStatistics", "--omitLocusTable", "--omitPerSampleStats", "--printBaseCounts", "o"=pathnameDT, verbose=verbose)

      verbose && cat(verbose, "GATK system result: ", res)
      verbose && exit(verbose)


      # (d) Cleanup GATK output file and resave
      verbose && enter(verbose, "Parsing and cleaning up the GATK output file")
      db <- TabularTextFile(pathnameDT)
      verbose && print(verbose, db)

      verbose && enter(verbose, "Reading")
      colClassPatterns <- c("(Locus|_base_counts)"="character")
      data <- readDataFrame(db, colClassPatterns=colClassPatterns)
      verbose && str(verbose, data)
      verbose && exit(verbose)

      # Drop SNPs with no coverage?
      if (dropEmpty) {
        verbose && enter(verbose, "Dropping SNPs with zero coverage")
        pattern <- "A:0 C:0 G:0 T:0 N:0"
        empty <- grep(pattern, data[[ncol(data)]], fixed=TRUE)
        verbose && printf(verbose, "Number of empty SNPs: %d (%g%%) of %d\n", length(empty), 100*length(empty)/nrow(data), nrow(data))
        if (length(empty) > 0L) {
          data <- data[-empty,]
        }
        rm(empty)
        verbose && str(verbose, data)
        verbose && exit(verbose)
      } #if (dropEmpty)

      verbose && enter(verbose, "Parsing genomic locations")
      # Parse (chromosome, position)
      chr <- gsub(":.*", "", data$Locus)
      pos <- gsub(".*:", "", data$Locus)
      pos <- as.integer(pos)
      .stop_if_not(length(pos) == length(chr))
      verbose && exit(verbose)

      verbose && enter(verbose, "Parsing (A,C,G,T,N) counts")
      # Parse (A,C,G,T) counts
      # NOTE: Some allele have non-zero "N" counts, which happens
      # when a read is aligned to a SNP, but its nucleotide at the
      # SNP position could not be called.  Because of this, we need
      # to keep the "N" column as well.
      counts <- data[[ncol(data)]]
      rm(data)
      counts <- gsub("[ACGTN]:", ":", counts)
      counts <- gsub(" ", "", counts, fixed=TRUE)
      counts <- gsub("^:", "", counts)
      counts <- strsplit(counts, split=":", fixed=TRUE)
      ns <- sapply(counts, FUN=length)
      .stop_if_not(all(ns == 5L))
      counts <- unlist(counts, use.names=FALSE)
      counts <- as.integer(counts)
      counts <- matrix(counts, ncol=5L, byrow=TRUE)
      colnames(counts) <- c("A", "C", "G", "T", "N")
      .stop_if_not(nrow(counts) == length(chr))
      verbose && exit(verbose)

      verbose && enter(verbose, "Summaries of nucleotide counts")
      countsT <- counts[,c("A", "C", "G", "T"),drop=FALSE]

      # Per-nucleotide coverage
      alleleCoverage <- colSums(countsT)
      alleleCoverage <- as.integer(alleleCoverage)
      names(alleleCoverage) <- colnames(countsT)

      # Total coverage
      totalCoverage <- sum(alleleCoverage)
      verbose && printf(verbose, "Total number of reads (with called alleles) covering a SNP: %d\n", totalCoverage)

      # Sanity check
      if (dropEmpty) {
        .stop_if_not(all(totalCoverage > 0L))
      }

      # Per-SNP coverage
      coverage <- rowSums(countsT, na.rm=TRUE)
      coverage <- as.integer(coverage)
      tblCoverage <- table(coverage)
      verbose && cat(verbose, "Distribution of SNP coverages:")
      verbose && print(verbose, tblCoverage)

      # Identify homozygous and heterozygous SNPs (with coverage >= 2)
      isHom <- rowAnys(counts == coverage)
      isHom[(coverage < 2L)] <- NA # Unknown
      isHet <- !isHom
      nbrOfHoms <- sum(isHom, na.rm=TRUE)
      nbrOfHets <- sum(isHet, na.rm=TRUE)
      verbose && printf(verbose, "Number of homozygous SNPs (with coverage >= 2): %d\n", nbrOfHoms)
      verbose && printf(verbose, "Number of heterozygous SNPs (with coverage >= 2): %d\n", nbrOfHets)

      tblHetCoverage <- table(isHet)
      verbose && cat(verbose, "Distribution of heterozygous SNP coverages:")
      verbose && print(verbose, tblHetCoverage)

      verbose && exit(verbose)


      verbose && enter(verbose, "Writing pruned GATK result file")
      header <- list(
        description = "GATK DepthOfCoverage Results",
        totalCoverage = totalCoverage,
        alleleCoverage = sprintf("%s:%d", names(alleleCoverage), alleleCoverage),
        nbrOfHoms = nbrOfHoms,
        nbrOfHets = nbrOfHets,
        tblCoverage = sprintf("%s:%d", names(tblCoverage), tblCoverage),
        tblHetCoverage = sprintf("%s:%d", names(tblHetCoverage), tblHetCoverage),
        dropEmpty = dropEmpty
      )

      # Updated allele count data
      data <- data.frame(chromosome=chr, position=pos, counts)
      rm(ns, chr, pos, counts)

      # Store
      pathnameDTT <- sprintf("%s.tmp", pathnameDT) # AD HOC
      pathnameDTT <- Arguments$getWritablePathname(pathnameDTT, mustNotExist=TRUE)

      pkg <- aroma.seq
      createdBy <- sprintf("%s v%s (%s)", getName(pkg), getVersion(pkg), getDate(pkg))
      writeDataFrame(data, file=pathnameDTT, header=header, createdBy=createdBy)
      rm(data, header)
      file.remove(pathnameDT) # AD HOC
      pathnameDT <- popTemporaryFile(pathnameDTT)
      verbose && exit(verbose)

      pathnameD <- popTemporaryFile(pathnameDT)
    } # if (overwrite || !isFile(...))

    # Parse GATK output file
    db <- TabularTextFile(pathnameD)
    verbose && print(verbose, db)

    verbose && exit(verbose)
  } # for (ii ...)

  res <- getOutputDataSet(this)
  verbose && print(verbose, res)

  verbose && exit(verbose)

  res
}) # process()



setMethodS3("readGatkCountFile", "GatkAlleleCounting", function(this, array, ..., verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose)
  if (verbose) {
    pushState(verbose)
    on.exit(popState(verbose))
  }

  verbose && enter(verbose, "Reading GATK count file")

  ds <- getOutputDataSet(this)
  array <- Arguments$getIndex(array, max=length(ds))

  df <- ds[[array]]
  verbose && print(verbose, df)

  # Parse GATK results
  pathname <- getPathname(df)
  db <- TabularTextFile(pathname)
  verbose && print(verbose, db)

  colClassPatterns <- c("*"="integer", "chromosome"="character")
  data <- readDataFrame(db, colClassPatterns=colClassPatterns)

  verbose && exit(verbose)

  data
}, protected=TRUE) # readGatkCountFile()


setMethodS3("getCombineBy", "GatkAlleleCounting", function(this, ...) {
  combineBy <- function(dfList) {
    # Drop loci with zero coverage
    countKeys <- c("A", "C", "G", "T", "N")
    dfList <- lapply(dfList, FUN=function(df) {
      dfT <- df[,countKeys]
      counts <- Reduce("+", dfT)
      df[counts > 0L,]
    })

    # Stack
    df <- Reduce(rbind, dfList)
    rm(dfList)

    dfK <- df[,c("chromosome", "position")]
    counts <- as.matrix(df[,countKeys])
    keys <- Reduce(function(...) paste(..., sep="\t"), dfK)
    dups <- duplicated(keys)

    # Unique (chromosome,position)
    dfU <- dfK[!dups,]
    keysU <- keys[!dups]
    countsU <- counts[!dups,,drop=FALSE]

    # Remainders
    keys <- keys[dups]
    counts <- counts[dups,,drop=FALSE]

    while(length(keys) > 0L) {
      dups <- duplicated(keys)

      # Unique (chromosome,position)
      keysT <- keys[!dups]
      countsT <- counts[!dups,,drop=FALSE]

      # Add to total counts
      idxs <- match(keysT, table=keys)
      .stop_if_not(all(is.finite(idxs)))
      countsU[idxs,] <- countsU[idxs,] + countsT

      # Remainders
      keys <- keys[dups]
      counts <- counts[dups,,drop=FALSE]
    } # while()

    data <- cbind(dfU, countsU)
    data
  } # combineBy()

  combineBy
}, protected=TRUE, static=TRUE)
