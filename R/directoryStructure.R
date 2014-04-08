# Pathnames:
# <rootpath>/<dataset>/<chiptype>/
#   pattern <- "(.*)/(.*)/(.*)/"
#   replacement <- c(rootpath="\\1", dataset="\\2", chiptype="\\3")
# <rootpath>/<dataset>/<organism>/*
#   pattern <- "(.*)/(.*)/(.*)/"
#   replacement <- c(rootpath="\\1", dataset="\\2", organism="\\3")
# <rootpath>/<dataset>/<organism>/<sample>/.*
#   pattern <- "(.*)/(.*)/(.*)/(.*)/"
#   replacement <- c(rootpath="\\1", dataset="\\2", organism="\\3", sample="\\4")


setMethodS3("directoryStructure", "NULL", function(struct, ...) NULL)

setMethodS3("directoryStructure", "list", function(struct, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'struct':
  names <- names(struct);
  for (name in c("pattern", "replacement")) {
    if (!is.element(name, names)) {
      throw(sprintf("List argument 'struct' does not have element %s: %s", sQuote(name), paste(sQuote(names), collapse=", ")));
    }
  }
  struct;
}) # directoryStructure()


setMethodS3("directoryStructure", "character", function(struct, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'struct':
  if (missing(struct)) {
    struct <- match.arg(struct);
  }
  struct <- Arguments$getCharacter(struct, length=c(1L, 1L));

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Create a directory structures?
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  parts <- strsplit(struct, split="/", fixed=TRUE);
  parts <- unlist(parts, use.names=FALSE);
  names <- gsub("^<(.*)>$", "\\1", parts);
  pattern <- paste(rep("([^/]*)", length=length(parts)), collapse="/");
  replacement <- sprintf("\\%d", seq_along(names));
  names(replacement) <- names;
  struct <- list(
    pattern = pattern,
    replacement = replacement
  );

  # Validate
  struct <- directoryStructure(struct, ...);

  struct;
}) # directoryStructure()

setMethodS3(".findDefaultDirectoryStructure", "GenericDataFile", function(this, ...) {
  fcn <- NULL;
  for (class in class(this)) {
    fcn <- getS3method("directoryStructure", class, optional=TRUE);
    if (!is.null(fcn)) {
      args <- formals(fcn);
      if (is.element("default", names(args))) return(args$default);
    }
  }
  throw(sprintf("Failed to locate default directory structure for class '%s'", class(this)[1L]));
})

setMethodS3(".findDefaultDirectoryStructure", "GenericDataFileSet", function(this, ...) {
  # Infer 'default' from corresponding file class.
  className <- this$getFileClass();
  clazz <- Class$forName(className);
  classNames <- class(newInstance(clazz));
  fcn <- NULL;
  for (class in classNames) {
    fcn <- getS3method("directoryStructure", class, optional=TRUE);
    if (!is.null(fcn)) {
      args <- formals(fcn);
      if (is.element("default", names(args))) return(args$default);
    }
  }
  throw(sprintf("Failed to locate default directory structure for class '%s'", className));
})

setMethodS3("directoryStructure", "GenericDataFile", function(this, default=NULL, ...) {
  parts <- this$.directoryStructure;
  if (is.null(parts)) {
    if (is.null(default)) default <- .findDefaultDirectoryStructure(this);
    parts <- directoryStructure(default, ...);
  }
  parts;
})

setMethodS3("directoryStructure", "GenericDataFileSet", function(this, default=NULL, ...) {
  parts <- this$.directoryStructure;
  if (is.null(parts)) {
    if (is.null(default)) default <- .findDefaultDirectoryStructure(this);
    parts <- directoryStructure(default, ...);
  }
  parts;
})

setMethodS3("directoryStructure<-", "GenericDataFile", function(this, ..., value) {
  if (missing(value)) { args <- list(...); value <- args[[length(args)]] };
  struct <- directoryStructure(value);
  this$.directoryStructure <- struct;
  invisible(this);
})

setMethodS3("directoryStructure<-", "GenericDataFileSet", function(this, ..., value) {
  if (missing(value)) { args <- list(...); value <- args[[length(args)]] };
  struct <- directoryStructure(value);
  this$.directoryStructure <- struct;
  # Update all files accordingly
  this <- updateDirectoryStructure(this);
  invisible(this);
})

setMethodS3("updateDirectoryStructure", "GenericDataFileSet", function(this, ...) {
  struct <- this$.directoryStructure;
  files <- this$files;
  files <- lapply(files, FUN=`directoryStructure<-`, struct);
  this$files <- files;
  invisible(this);
}, protected=TRUE)


setMethodS3("directoryItems", "character", function(paths, struct, ..., as="list") {
  # Argument 'struct':
  struct <- directoryStructure(struct, ...);

  # Nothing to do?
  if (is.null(struct)) return(list());

  # Append optional slush/tail to the end
  pattern <- sprintf("%s(|/(.*))", struct$pattern);
  tail <- sprintf("\\%d", length(struct$replacement)+1L);
  names(tail) <- "<tail>";
  replacement <- c(struct$replacement, tail);

  # Parse path according to 'struct'
  paths <- gsub("\\", "/", paths, fixed=TRUE);
  res <- agsub(pattern=pattern, replacement=replacement, paths, ..., as=as);

  res;
}) # directoryItems()


setMethodS3("directoryItems", "GenericDataFile", function(this, ...) {
  struct <- directoryStructure(this);
  # Nothing to do?
  if (is.null(struct)) return(list());

  pathname <- getPathname(this);
  directoryItems(pathname, struct=struct, ...);
}, protected=TRUE)


setMethodS3("directoryItems", "GenericDataFileSet", function(this, ...) {
  struct <- directoryStructure(this);
  # Nothing to do?
  if (is.null(struct)) return(list());

  path <- getPath(this);
  pathname <- file.path(path, NA_character_);
  directoryItems(pathname, struct=struct, ...);
}, protected=TRUE)


setMethodS3("directoryItem", "GenericDataFileSet", function(this, name, default=NULL, ..., mustExist=TRUE) {
  items <- directoryItems(this, ...);
  if (!is.element(name, names(items))) {
    if (is.null(default)) {
      if (mustExist) {
        path <- getPath(this);
        throw(sprintf("Cannot infer %s from path %s using set directory structure (%s).", sQuote(name), sQuote(path), paste(sQuote(names(items)), collapse=", ")));
      }
    }
    items[[name]] <- default;
  }
  value <- items[[name]];
  value;
}, protected=TRUE) # directoryItem()

setMethodS3("directoryItem", "GenericDataFile", function(this, name, default=NULL, ..., mustExist=TRUE) {
  items <- directoryItems(this, ...);
  if (!is.element(name, names(items))) {
    if (is.null(default)) {
      if (mustExist) {
        pathname <- getPathname(this);
        throw(sprintf("Cannot infer %s from pathname %s using set directory structure (%s).", sQuote(name), sQuote(pathname), paste(sQuote(names(items)), collapse=", ")));
      }
    }
    items[[name]] <- default;
  }
  value <- items[[name]];
  value;
}, protected=TRUE) # directoryItem()


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AD HOC
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("getOrganismName", "GenericDataFile", function(this, ...) {
  directoryItem(this, name="organism");
})

setMethodS3("getDataSetName", "GenericDataFile", function(this, ...) {
  directoryItem(this, name="dataset");
})

setMethodS3("getSampleName", "GenericDataFile", function(this, ...) {
  directoryItem(this, name="sample", default=getFullName(this, ...));
})


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# AROMA.SEQ GENERIC
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
setMethodS3("directoryStructure", "AromaSeqDataFile", function(this, default="<rootpath>/<dataset>/<organism>/<sample>/", ...) {
  if (is.null(default)) default <- .findDefaultDirectoryStructure(this);
  NextMethod("directoryStructure", default=default);
})

setMethodS3("directoryStructure", "AromaSeqDataFileSet", function(this, default="<rootpath>/<dataset>/<organism>/<sample>/", ...) {
  if (is.null(default)) default <- .findDefaultDirectoryStructure(this);
  NextMethod("directoryStructure", default=default);
})

setMethodS3("getOrganism", "AromaSeqDataFile", function(this, ...) {
  directoryItem(this, name="organism");
})

setMethodS3("getOrganism", "AromaSeqDataFileSet", function(this, ...) {
  directoryItem(this, name="organism");
})

setMethodS3("getDefaultFullName", "AromaSeqDataFile", function(this, ...) {
  value <- directoryItem(this, name="sample", mustExist=FALSE);
  if (is.null(value)) {
    value <- NextMethod("getDefaultFullName");
  } else {
    pattern <- getExtensionPattern(this);
    value <- gsub(pattern, "", value);
  }
  value;
})

setMethodS3("getDefaultFullName", "AromaSeqDataFileSet", function(this, ...) {
  directoryItem(this, name="dataset");
})


############################################################################
# HISTORY
# 2014-04-07
# o ROBUSTNESS: Now directory structure items may not contain slashes.
# o Added argument 'firstOnly' to directoryItems() for GenericDataFileSet.
# o CLEANUP: Added .findDefaultDirectoryStructure().
# 2013-11-10
# o Added directoryStructure() for BAM, SAM and FASTQ classes.
# o Added directoryStructure() and ditto replacement functions,
#   directoryItem()/directoryItems() for character, GenericDataFile
#   and GenericDataFileSet.
# o Created.
############################################################################
