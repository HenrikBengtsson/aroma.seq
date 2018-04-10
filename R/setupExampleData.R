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
setupExampleData <- function(dirs=c("annotationData", "bamData", "fastqData*"), ...) {
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

    # Nothing to do?
    if (!isDirectory(dirS)) next;

    if (link) {
      dir <- Arguments$getWritablePath(dir);
      for (dirSS in list.files(path=dirS)) {
        createLink(link=file.path(dir, dirSS), target=file.path(dirS, dirSS));
      }
    } else {
      copyDirectory(dirS, to=dir, copy.mode=FALSE, overwrite=FALSE);
    }

    # Sanity check
    stop_if_not(isDirectory(dir));
  } # for (ii ...)


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Data set: GATK Resource Bundle (iff available)
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (isCapableOf(aroma.seq, "gatk")) {
    bin <- findGATK();
    srcPath <- file.path(dirname(bin), "resources");
    srcPath <- Arguments$getReadablePath(srcPath, mustExist=FALSE);
    if (isDirectory(srcPath)) {
      dataset <- "GATKResourceBundle";
      organism <- "GATKExample";

      if (is.element("annotationData", dirs)) {
        path <- file.path("annotationData", "organisms", organism);
        path <- Arguments$getWritablePath(path);
        pathnames <- dir(path=srcPath, pattern="^exampleFASTA[.]", full.names=TRUE);
        sapply(pathnames, FUN=function(pathname) {
          copyFile(pathname, file.path(path, basename(pathname)), copy.mode=FALSE, skip=TRUE);
        })
      }

      if (is.element("bamData", dirs)) {
        path <- file.path("bamData", dataset, organism);
        path <- Arguments$getWritablePath(path);

        pathnames <- dir(path=srcPath, pattern="^exampleBAM[.]", full.names=TRUE)
        sapply(pathnames, FUN=function(pathname) {
          copyFile(pathname, file.path(path, basename(pathname)), copy.mode=FALSE, skip=TRUE)
        })
      }
    } # if (isDirectory(srcPath))
  }

  invisible(dirs);
} # setupExampleData()
