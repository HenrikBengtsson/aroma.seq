###########################################################################/**
# @RdocClass BamDataFile
#
# @title "The abstract BamDataFile class"
#
# \description{
#  @classhierarchy
#
#  A BamDataFile object represents a BAM file.
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
# @author
#
# \references{
#  [1] The SAM Format Specification Working Group,
#      \emph{The SAM Format Specification}, Sept 7, 2011.\cr
# }
#
# \seealso{
#   An object of this class is typically part of an 
#   @see "BamDataSet".
# }
#*/###########################################################################
setConstructorS3("BamDataFile", function(...) {
  extend(GenericDataFile(...), "BamDataFile");
})


setMethodS3("as.character", "BamDataFile", function(x, ...) {
  # To please R CMD check
  this <- x;

  s <- NextMethod("as.character", ...);
  class <- class(s);

  s <- c(s, sprintf("Has index file (*.bai): %s", hasIndex(this)));
  s <- c(s, sprintf("Is sorted: %s", isSorted(this)));
  if (hasIndex(this)) {
    n <- nbrOfTargets(this);
    s <- c(s, sprintf("Number of targets: %s", n));
    len <- getTotalTargetLength(this);
    s <- c(s, sprintf("Total target length: %.3gMb (%.0f bases)", len/1e9, len));
    names <- getTargetNames(this);
    s <- c(s, sprintf("Targets: [%d] %s", n, hpaste(names)));
  }

  class(s) <- class;
  s;
})

setMethodS3("buildIndex", "BamDataFile", function(this, ..., skip=!overwrite, overwrite=FALSE) {
  pathname <- getPathname(this);
  pathnameBAI <- sprintf("%s.bai", pathname);

  if (hasIndex(this)) {
    if (skip) {
      return(invisible(pathnameBAI));
    }
    if (!overwrite) {
      throw("Cannot build index file (*.bai). File already exists: ", pathnameBAI);
    }
  }

  pathnameT <- Rsamtools::indexBam(pathname);

  invisible(pathnameBAI);
})


setMethodS3("hasIndex", "BamDataFile", function(this, ...) {
  pathname <- getPathname(this);
  pathnameBAI <- sprintf("%s.bai", pathname);
  isFile(pathnameBAI);
})

# \details{
#   BAM headers typically contain an \code{"@HD VN:1.0 SO:<value>"} entry,
#   where \code{<value>} indicates whether the aligned reads are sorted 
#   or not.  Unfortunately, this entry is neither enforced nor has it to
#   be correct [1,2].
#
#   Instead, we consider a BAM file to be sorted if and only if it has
#   an index file.  The rationale is that it is not possible to index
#   a BAM file unless it is sorted first.
# }
#
# \references{
#   [1] Question: is my BAM file sorted?, Biostar, 2011,
#   \url{http://www.biostars.org/post/show/5256/is-my-bam-file-sorted/}\cr
#   [2] Asking for suggestiona on samtools bug fixing, SEQanswers, 2010,
#   \url{http://seqanswers.com/forums/showthread.php?t=3739}\cr
# }
setMethodS3("isSorted", "BamDataFile", function(this, ...) {
  isTRUE(hasIndex(this));
})


# Argument '...' must be 2nd to match the generic base::sort() function.
setMethodS3("sort", "BamDataFile", function(x, ..., force=FALSE) {
  # To please R CMD check
  this <- x;

  # Nothing todo?
  if (!force && isSorted(this)) {
    return(this);
  }

  throw("Not yet implemented!");
})


setMethodS3("nbrOfSeqs", "BamDataFile", function(this, ...) {
  nbrOfTargets(this);
})

setMethodS3("getTargets", "BamDataFile", function(this, ...) {
  hdr <- getHeader(this);
  targets <- hdr$targets;
  targets;  
})

setMethodS3("nbrOfTargets", "BamDataFile", function(this, ...) {
  length(getTargets(this));
})

setMethodS3("getTargetNames", "BamDataFile", function(this, ...) {
  names(getTargets(this));
})

setMethodS3("getTargetLengths", "BamDataFile", function(this, ...) {
  getTargets(this);
})

setMethodS3("getTotalTargetLength", "BamDataFile", function(this, ...) {
  sum(as.numeric(getTargets(this)));
})

setMethodS3("getHeader", "BamDataFile", function(this, force=FALSE, ...) {
  header <- this$.header;
  if (force || is.null(header)) {
    header <- readHeader(this, ...);
    this$.header <- header;
  }
  header;
})

setMethodS3("readHeader", "BamDataFile", function(this, ...) {
  pathname <- getPathname(this);
  if (!hasIndex(this)) {
    throw("Cannot read header. Index file (*.bai) is missing: ", pathname);
  }
  bf <- Rsamtools::BamFile(pathname);
  hdr <- Rsamtools::scanBamHeader(bf);
  hdr;
}, private=TRUE)


# \section{Read Groups (RG) specification}{
#   The RG fields populated are ID, SM, PL and LB, where
#   ID=identifier, SM=sample, PL=platform/technology, and
#   LB=library (DNA library prep identifier).
#   ID: Each RG must have a unique ID.
#   SM: Sample. Use pool name where a pool is being sequenced.
#   PM: Platform/technology used to produce the reads. Valid values: 
#       CAPILLARY, HELICOS, ILLUMINA, IONTORRENT, LS454, PACBIO, and SOLID.
#   LB: Library.
#
#   Note that tools such as GATK (Broad Institute) requires RG:s.
# }
#
# \references{
#  [1] The SAM Format Specification Working Group,
#      \emph{The SAM Format Specification}, Sept 7, 2011.\cr
#  [2] "What is '@RG ID' versus '@RG SM'"
#      \url{http://seqanswers.com/forums/showthread.php?t=9784}\cr
# }
setMethodS3("getReadGroups", "BamDataFile", function(this, ...) {
  hdr <- readHeader(this, ...);
  text <- hdr$text;
  keys <- names(text);
  rgList <- text[(keys == "@RG")];
  rgList <- lapply(rgList, FUN=function(params) {
    pattern <- "([^:]*):(.*)";
    keys <- gsub(pattern, "\\1", params);
    values <- gsub(pattern, "\\2", params);
    names(values) <- keys;
    values;
  });
  rgList;
})


# @RdocMethod replaceAllReadGroups
# @title "Writes a new BAM file with all existing read groups replaced by one new read group"
#
# \description{
#   @get "title".
# }
#
# \arguments{
#  \item{sample}{Specifies the \code{SM} read group.}
#  \item{library}{Specifies the \code{LB} read group.}
#  \item{platform}{Specifies the \code{PL} read group.}
#  \item{platformUnit}{Specifies the \code{PU} read group.}
# }
setMethodS3("replaceAllReadGroups", "BamDataFile", function(this, sample="*", library="*", platform="*", platformUnit="*", sequencingCenter="*", description="*", runDate="*", ..., validate=TRUE, skip=!overwrite, overwrite=FALSE, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Local functions
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  getReadGroupsOnce <- function(rg, ...) {
    if (is.null(rg)) {
      rgList <- getReadGroups(this);
      names(rgList) <- NULL;
      rg <- unlist(rgList);
      rg <- rev(rg);
      dups <- duplicated(names(rg));
      rg <- rg[!dups];
      rg <- rev(rg);
    }
    rg;
  } # getReadGroupsOnce()



  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'sample':
  sample <- Arguments$getCharacter(sample);

  # Argument 'library':
  library <- Arguments$getCharacter(library);

  # Argument 'platform':
  platform <- Arguments$getCharacter(platform);

  # Argument 'platformUnit':
  platformUnit <- Arguments$getCharacter(platformUnit);

  # Argument 'sequencingCenter':
  sequencingCenter <- Arguments$getCharacter(sequencingCenter);

  # Argument 'description':
  description <- Arguments$getCharacter(description);

  # Argument 'runDate':
  if (inherits(runDate, "Date")) {
    runDate <- Arguments$getCharacter(runDate);
  } else if (inherits(runDate, "POSIXct") || inherits(runDate, "POSIXlt")) {
    runDate <- Arguments$getCharacter(runDate);
  } else {
    runDate <- Arguments$getCharacter(runDate);
  }

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }

  verbose && enter(verbose, "Updating Read Groups");
 

  pathname <- getPathname(this);

  # Output pathname
  pathnameD <- gsub("[.]bam$", ",RG.bam", pathname);
  pathnameD <- Arguments$getWritablePathname(pathnameD, mustNotExist=(!skip && !overwrite));

  if (isFile(pathnameD)) {
    if (skip) {
      verbose && cat(verbose, "Already exists. Skipping.");
      bf <- newInstance(this, pathnameD);
      buildIndex(bf, skip=TRUE, verbose=less(verbose, 10));
      verbose && exit(verbose);
      return(bf);
    }
    if (!overwrite) {
      throw("Cannot replace SAM read groups. File already exists: ", pathnameD);
    }
  }


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  optional <- c(CN="sequencingCenter", DS="description", DT="runDate");
  map <- c(SM="sample", LB="library", PL="platform", PU="platformUnit",
           optional);
  rgList <- list();
  for (key in names(map)) {
    varname <- map[key];
    rgList[[key]] <- get(varname, inherits=FALSE);
  }

  # Append user arguments '...'
  args <- list(...);
  for (key in names(args)) {
    rgList[[key]] <- args[[key]];
  }
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, rgList);

  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, rgList);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Use default values from existing read groups of the input file?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  rg <- NULL;
  for (key in names(rgList)) {
    value <- rgList[[key]];

    # Use default?
    if (value == "*") {
      rg <- getReadGroupsOnce(rg);
      value <- rg[key];
      if (is.na(value)) {
        if (!is.element(key, names(optional))) {
          throw(sprintf("Read group '%s' could not be inferred from input file ('%s'). The following read groups was found: %s", key, pathname, hpaste(rg)));
        }
        value <- NULL;
      }
      rgList[[key]] <- value;
    }
  } # for (key ...)
  verbose && cat(verbose, "Arguments:");
  verbose && str(verbose, rgList);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate read groups?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (validate) {
    # Validate platform
    knownPlatforms <- c("CAPILLARY", "HELICOS", "ILLUMINA", "IONTORRENT", "LS454", "PACBIO", "SOLID");
    platform <- rgList[["PL"]];
    if (!is.element(toupper(platform), knownPlatforms)) {
      throw("Unknown 'PL' (platform) value: ", platform);
    }

    # Validate run date
    runDate <- rgList[["DT"]];
    # TODO...
  } # if (validate)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Write new BAM file
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  args <- list("AddOrReplaceReadGroups", I=pathname, O=pathnameD);
  names(rgList) <- sprintf("RG%s", names(rgList));
  args <- c(args, rgList);
  args$verbose <- less(verbose, 10);
  verbose && str(verbose, args);

  res <- do.call(systemPicard, args);

  bf <- newInstance(this, pathnameD);
  buildIndex(bf, overwrite=TRUE, verbose=less(verbose, 10));

  verbose && exit(verbose);

  bf;
}) # replaceAllReadGroups()


############################################################################
# HISTORY:
# 2012-06-28
# o Added getReadGroups() and replaceAllReadGroups() and buildIndex().
# o Created.
############################################################################
