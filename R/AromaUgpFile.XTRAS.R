setMethodS3("anyDuplicated", "AromaUgpFile", function(x, ...) {
  # To please R CMD check
  this <- x;
  data <- readDataFrame(this, ...);
  anyDuplicated(data, ...);
})


setMethodS3("writeBedDataFile", "AromaUgpFile", function(this, ..., path=getPath(this), chrMap=NULL, skip=TRUE, overwrite=!skip, verbose=FALSE) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  # Argument 'chrMap':
  if (!is.null(chrMap)) {
    chrMap <- Arguments$getIndices(chrMap, useNames=TRUE);
    if (is.null(names(chrMap))) {
      throw("Argument 'chrMap' should have names.");
    }
  }

  # Argument 'skip':
  skip <- Arguments$getLogical(skip);

  # Argument 'overwrite':
  overwrite <- Arguments$getLogical(overwrite);

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enter(verbose, "Exporting UGP as BED file");
  filename <- sprintf("%s.bed", getFilename(this));
  pathname <- Arguments$getReadablePathname(filename, path=path, mustExist=FALSE);
  verbose && cat(verbose, "UGP pathname: ", getPathname(this));
  verbose && cat(verbose, "BED pathname: ", pathname);

  # Already done?
  if (skip && isFile(pathname)) {
    verbose && cat(verbose, "Already exported. Skipping.");
    res <- GenericDataFile(pathname);
    verbose && print(verbose, res);
    verbose && exit(verbose);
    return(invisible(res));
  }

  pathname <- Arguments$getWritablePathname(pathname, mustNotExist=!overwrite);

  verbose && enter(verbose, "Reading data");
  data <- readDataFrame(this, verbose=less(verbose, 10));
  verbose && str(verbose, data);
  verbose && exit(verbose);

  if (!is.null(chrMap)) {
    verbose && enter(verbose, "Renaming chromosomes");
    verbose && print(verbose, chrMap);
    keys <- names(chrMap);
    chrs <- data$chromosome;
    for (kk in seq(along=keys)) {
      key <- keys[kk];
      chr <- chrMap[kk];
      idxs <- which(chrs == chr);
      if (length(idxs) > 0L) {
        chrs[idxs] <- key;
      }
    }
    data$chromosome <- chrs;
    # Not needed anymore
    rm(keys, chrs);
    verbose && exit(verbose);
  }

  verbose && enter(verbose, "Writing BED file");
  pathnameD <- writeDataFrame(data, file=pathname, col.names=FALSE, header=NULL);
  verbose && exit(verbose);

  # Not needed anymore
  rm(data);

  res <- GenericDataFile(pathname);
  verbose && print(verbose, res);

  verbose && exit(verbose);

  invisible(res);
}) # writeBedDataFile()


############################################################################
# HISTORY:
# 2012-10-31.
# o Added anyDuplicated() for AromaUgpFile.
# o Added writeBedDataFile() for AromaUgpFile.
# o Created.
############################################################################
