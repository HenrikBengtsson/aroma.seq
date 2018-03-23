.findCache <- local({
  db <- new.env(parent=emptyenv());

  getKey <- function(name, version=NULL) {
    if (!is.null(version)) {
      version <- as.character(version);
      version <- rep(version, length.out=2L);
      name <- sprintf("%s;version=[%s,%s]", name, version[1L], version[2L]);
    }
    name;
  } # getKey()

  function(name, version=NULL, path, ...) {
    if (missing(path)) {
      if (missing(name)) {
        return(as.list(db));
      }
      key <- getKey(name, version=version);
      res <- db[[key]];
      if (!is.null(res) || !is.null(version)) {
        return(res);
      }
      # Find all possible matches
      pattern <- sprintf("^%s;version=", key);
      res <- as.list(db);
      keys <- grep(pattern, names(res));
      res <- res[keys];
      if (length(res) == 0L) res <- NULL;
      return(res);
    } else {
      key <- getKey(name, version=version);
      db[[key]] <- list(key=key, name=name, version=version, path=path);
    }
    invisible(NULL);
  }
}) # .findCache()
