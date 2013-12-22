setMethodS3("loadCacheFile", "GenericDataFile", function(this, removeOldCache=TRUE, ...) {
  # Local functions
  cleanup <- function() {
    if (removeOldCache && isFile(pathname)) file.remove(pathname);
  }

  pathnameS <- getPathname(this);
  pathname <- sprintf("%s.Rcache", pathnameS);
  cache <- loadCache(pathname=pathname);

  # Not the same class?
  if (!identical(cache$class, class(this)[1L])) {
    cleanup();
    return(NULL);
  }

  # No checksum?
  checksum <- cache$checksum;
  if (is.null(checksum)) {
    cleanup();
    return(NULL);
  }

  # Non-match checksum?
  cf <- getChecksumFile(this);
  checksumF <- readChecksum(cf);
  if (!identical(checksum, checksumF)) {
    cleanup();
    return(NULL);
  }

  cache;
}, protected=TRUE)


setMethodS3("loadCacheFileItem", "GenericDataFile", function(this, name, ...) {
  cache <- loadCacheFile(this, ...);
  cache[[name]];
}, protected=TRUE)



setMethodS3("saveCacheFile", "GenericDataFile", function(this, cache=NULL, ...) {
  # Update/set class
  cache$class <- class(this)[1L];

  # Update/set/generate checksum
  cf <- getChecksumFile(this);
  checksum <- readChecksum(cf);
  cache$checksum <- checksum;

  pathnameS <- getPathname(this);
  pathname <- sprintf("%s.Rcache", pathnameS);
  saveCache(cache, pathname=pathname);
}, protected=TRUE)


setMethodS3("saveCacheFileItem", "GenericDataFile", function(this, name=NULL, value=NULL, ...) {
  args <- list();
  if (!is.null(name)) {
    args[[name]] <- value;
  }
  args <- c(args, list(...));
  nargs <- length(args);

  # Nothing to do?
  if (nargs == 0L) {
    return();
  }
  names <- names(args);
  if (is.null(names)) {
    throw("Arguments to saveCacheFileItem() must be named.");
  }

  # Update existing/add to cache
  cache <- loadCacheFile(this);
  if (is.null(cache)) cache <- list();

  for (name in names) {
    cache[[name]] <- args[[name]];
  }

  saveCacheFile(this, cache=cache);
}, protected=TRUE)


setMethodS3("memoizedCall2", "GenericDataFile", function(this, what, ..., envir=parent.frame(), force=FALSE, verbose=FALSE) {
  # 1. Look for memoized results
  key <- list(method=what, ...);
  keyCS <- getChecksum(key);
  if (!force) {
    res <- loadCacheFileItem(this, name=keyCS);
    if (!is.null(res)) {
      if (verbose) cat("Returning cached results!");
      value <- res$value;
      return(value);
    }
  }

  # 2. Otherwise, call method with arguments
  args <- list(this, ...);
  value <- do.call(what, args=args, quote=FALSE, envir=envir);

  # 3. Memoize results
  saveCacheFileItem(this, name=keyCS, value=list(key=key, value=value));

  # 4. Return results
  value;
}, protected=TRUE)


############################################################################
# HISTORY:
# 2012-12-21
# o Created.
############################################################################
