.findCache <- local({
  db <- new.env(parent=emptyenv());

  getKey <- function(name, version=NULL) {
    if (!is.null(version)) {
      name <- sprintf("%s_v%s", name, as.character(version));
    }
    name;
  } # getKey()

  function(name, version=NULL, path, ...) {
    if (missing(path)) {
      if (missing(name)) {
        return(as.list(db));
      }
      key <- getKey(name, version=version);
      return(db[[key]]);
    } else {
      key <- getKey(name, version=version);
      db[[key]] <- list(key=key, name=name, version=version, path=path);
    }
    invisible(NULL);
  }
}) # .findCache()


############################################################################
# HISTORY:
# 2013-04-01
# o Created.
############################################################################
