setConstructorS3("FileGroupsInterface", function(...) {
  extend(Interface(), "FileGroupsInterface");
})


setMethodS3("getGroupBy", "FileGroupsInterface", function(this, ...) {
  params <- getParameters(this);
  params$groupBy;
}, protected=TRUE)


setMethodS3("getInputDataSet", "FileGroupsInterface", function(...) {
  NextMethod("getInputDataSet");
}, protected=TRUE)


setMethodS3("getGroups", "FileGroupsInterface", function(this, ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Identify groups
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  ds <- getInputDataSet(this);
  fullnames <- getFullNames(ds);
  # FIXME: Dropping "R1" suffix should really be done by the data set/files.
  fullnames <- sub("_(1|R1)$", "", fullnames);

  groups <- getGroupBy(this);
  if (is.null(groups)) {
    groups <- as.list(seq_along(ds));
    names(groups) <- fullnames;
  } else if (is.character(groups)) {
    if (groups == "name") {
      names <- getNames(ds);
      namesU <- unique(names);
      groups <- lapply(namesU, FUN=function(name) which(names == name));
      names(groups) <- namesU;
    }
  }
  # Sanity check
  stopifnot(is.list(groups));


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Add names to groups and indices
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Group names
  if (is.null(names(groups))) {
    names(groups) <- sprintf("Group_%d", seq_along(groups));
  }

  # Add index names, iff missing
  groups <- lapply(groups, FUN=function(idxs) {
    if (is.null(names(idxs))) {
      names(idxs) <- fullnames[idxs];
    }
    idxs;
  })


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validation
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Range and uniqueness check
  max <- length(ds);
  for (gg in seq_along(groups)) {
    idxs <- groups[[gg]];
    idxs <- Arguments$getIndices(idxs, max=max);
    dups <- duplicated(idxs);
    if (any(dups)) {
      throw(sprintf("Detected duplicated file indices in group %s: %s", names(groups)[gg], hpaste(idxs[dups])));
    }
  } # for (gg ...)

  # Additional class-specific validation, iff any
  validateGroups(this, groups);

  groups;
}, protected=TRUE) # getGroups()


setMethodS3("nbrOfGroups", "FileGroupsInterface", function(this, ...) {
  length(getGroups(this));
}, protected=TRUE)

setMethodS3("getGroupNames", "FileGroupsInterface", function(this, ...) {
  groups <- getGroups(this);
  names(groups);
}, protected=TRUE)

setMethodS3("validateGroups", "FileGroupsInterface", function(this, ...) {
  this;
})


############################################################################
# HISTORY:
# 2014-01-16 [HB]
# o Using this Interface for BamMerger and TopHat2Alignment.
# o Added FileGroupsInterface.
# o Created.
# Previous related history below:
# 2014-01-16 [HB]
# o BUG FIX: getGroups() of TopHat2Alignment would not generate the
#   correct sets of indices.
# o ROBUSTNESS: Now getGroups() of TopHat2Alignment assert that the
#   file indices identified for each group/sample is unique.
# 2013-11-22
# o Implemented getGroups() for BamMerger.
############################################################################
