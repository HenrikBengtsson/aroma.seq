###########################################################################/**
# @RdocFunction setupExampleData
#
# @title "Setups example data in the current directory"
#
# \description{
#  @get "title".
# }
#
# @synopsis
#
# \arguments{
#   \item{dirs}{A @character @vector specifying which directories to setup.}
#   \item{...}{Not used.}
# }
#
# @author "HB"
#
# @keyword internal
#*/###########################################################################
setupExampleData <- function(dirs=c("annotationData", "fastqData*"), ...) {
  # Argument 'dirs':
  dirs <- Arguments$getCharacters(dirs);

  links <- (regexpr("[*]$", dirs) != -1L);
  dirs <- gsub("[*]$", "", dirs);

  pathS <- system.file("exData", package="aroma.seq", mustWork=TRUE);

  for (ii in seq_along(dirs)) {
    dir <- dirs[ii];
    link <- links[ii];

    # Nothing to do?
    if (isDirectory(dir)) next;

    # Create
    dirS <- file.path(pathS, dir);
    if (link) {
      dir <- Arguments$getWritablePath(dir);
      for (dirSS in list.files(path=dirS)) {
        createLink(link=file.path(dir, dirSS), target=file.path(dirS, dirSS));
      }
    } else {
      copyDirectory(dirS, to=dir, overwrite=FALSE);
    }

    # Sanity check
    stopifnot(isDirectory(dir));
  } # for (ii ...)

  invisible(dirs);
} # setupExampleData()


############################################################################
# HISTORY:
# 2013-11-01
# o Added setupExampleData().
# o Created.
############################################################################
