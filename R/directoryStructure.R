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
  pattern <- paste(rep("(.*)", length=length(parts)), collapse="/");
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


setMethodS3("directoryStructure", "GenericDataFile", function(this, default=NULL, ...) {
  parts <- this$.directoryStructure;
  if (is.null(parts)) parts <- directoryStructure(default, ...);
  parts;
})

setMethodS3("directoryStructure", "GenericDataFileSet", function(this, default=NULL, ...) {
  parts <- this$.directoryStructure;
  if (is.null(parts)) parts <- directoryStructure(default, ...);
  parts;
})

setMethodS3("directoryStructure<-", "GenericDataFile", function(this, value) {
  struct <- directoryStructure(value);
  this$.directoryStructure <- struct;
  invisible(this);
})

setMethodS3("directoryStructure<-", "GenericDataFileSet", function(this, value) {
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

  # Parse path according to 'struct'
  paths <- gsub("\\", "/", paths, fixed=TRUE);
  res <- agsub(pattern=struct$pattern, replacement=struct$replacement, paths, ..., as=as);

  res;
}) # directoryItems()


setMethodS3("directoryItems", "GenericDataFile", function(this, ...) {
  path <- getPath(this);
  struct <- directoryStructure(this);
  directoryItems(path, struct=struct, ...);
}, protected=TRUE)


setMethodS3("directoryItems", "GenericDataFileSet", function(this, ...) {
  paths <- sapply(this, FUN=getPath);
  struct <- directoryStructure(this);
  directoryItems(paths, struct=struct, ...);
}, protected=TRUE)


setMethodS3("directoryItem", "GenericDataFileSet", function(this, name, default=NULL, ...) {
  items <- directoryItems(this, ...);
  if (!is.element(name, names(items))) {
    if (is.null(default)) {
      path <- getPath(this);
      throw(sprintf("Cannot infer %s from path %s using set directory structure (%s).", sQuote(name), sQuote(path), paste(sQuote(names(items)), collapse=", ")));
    }
    items[[name]] <- default;
  }
  value <- items[[name]];
  value;
}, protected=TRUE) # directoryItem()

setMethodS3("directoryItem", "GenericDataFile", function(this, name, default=NULL, ...) {
  items <- directoryItems(this, ...);
  if (!is.element(name, names(items))) {
    if (is.null(default)) {
      path <- getPath(this);
      throw(sprintf("Cannot infer %s from path %s using set directory structure (%s).", sQuote(name), sQuote(path), paste(sQuote(names(items)), collapse=", ")));
    }
    items[[name]] <- default;
  }
  value <- items[[name]];
  value;
}, protected=TRUE) # directoryItem()


setMethodS3("getOrganismName", "GenericDataFile", function(this, ...) {
  directoryItem(this, name="organism");
})

setMethodS3("getDataSetName", "GenericDataFile", function(this, ...) {
  directoryItem(this, name="dataset");
})

setMethodS3("getSampleName", "GenericDataFile", function(this, ...) {
  directoryItem(this, name="sample", default=getFullName(this));
})



setMethodS3("directoryStructure", "FastqDataFile", function(this, default="<rootpath>/<dataset>/<organism>/", ...) {
  NextMethod("directoryStructure", default=default);
})

setMethodS3("directoryStructure", "FastqDataSet", function(this, default=NULL, ...) {
  if (is.null(default)) {
    # Infer 'default' from corresponding file class.
    className <- this$getFileClass();
    fcn <- getS3method("directoryStructure", className);
    default <- formals(fcn)$default;
  }
  NextMethod("directoryStructure", default=default);
})


setMethodS3("directoryStructure", "BamDataFile", function(this, default="<rootpath>/<dataset>/<organism>/", ...) {
  NextMethod("directoryStructure", default=default);
})

setMethodS3("directoryStructure", "BamDataSet", function(this, default=NULL, ...) {
  if (is.null(default)) {
    # Infer 'default' from corresponding file class.
    className <- this$getFileClass();
    fcn <- getS3method("directoryStructure", className);
    default <- formals(fcn)$default;
  }
  NextMethod("directoryStructure", default=default);
})


setMethodS3("directoryStructure", "SamDataFile", function(this, default="<rootpath>/<dataset>/<organism>/", ...) {
  NextMethod("directoryStructure", default=default);
})

setMethodS3("directoryStructure", "SamDataSet", function(this, default=NULL, ...) {
  if (is.null(default)) {
    # Infer 'default' from corresponding file class.
    className <- this$getFileClass();
    fcn <- getS3method("directoryStructure", className);
    default <- formals(fcn)$default;
  }
  NextMethod("directoryStructure", default=default);
})


############################################################################
# HISTORY
# 2013-11-10
# o Added directoryStructure() for BAM, SAM and FASTQ classes.
# o Added directoryStructure() and ditto replacement functions,
#   directoryItem()/directoryItems() for character, GenericDataFile
#   and GenericDataFileSet.
# o Created.
############################################################################
