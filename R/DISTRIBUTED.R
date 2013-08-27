
#
# @RdocFunction .usePackage
#
# @title "Attaches a package"
#
# \description{
#  @get "title".  If not installed, it will be installed from one of
#  the known repositories.
# }
#
# @synopsis
#
# \seealso{
#   @see "base::library" and "base::install.packages".
# }
.usePackage <- function(pkg, minVersion=NULL, how=c("attach", "import"), install=FALSE, ...) {
  # Argument 'how':
  how <- match.arg(how);

  if (!is.null(minVersion) && !is.na(minVersion)) {
    ver <- packageVersion(pkg);
    if (ver < minVersion) {
      if (install) {
        install.packages(pkg, ...);
        res <- .usePackage(pkg, minVersion=minVersion, how=how, install=FALSE, ...);
        return(invisible(res));
      } else {
        throw(sprintf("%s v%s or newer is not installed: %s", pkg, minVersion, ver));
      }
    }
  }

  if (how == "attach") {
    require(pkg, character.only=TRUE) || throw("Package not attached: ", pkg);
  } else if (how == "import") {
    requireNamespace(pkg) || throw("Package not loaded: ", pkg);
  }

  res <- packageVersion(pkg);
  invisible(res);
} # .usePackage()


.useBatchJobs <- function(minVersion="*", ...) {
  if (identical(minVersion, "*")) {
    desc <- packageDescription("aroma.seq");
    desc <- desc[c("Depends", "Imports", "Suggests")];
    desc <- unlist(desc, use.names=FALSE);
    desc <- strsplit(desc, split=",", fixed=TRUE);
    desc <- unlist(desc, use.names=FALSE);
    desc <- gsub("\n", "", desc, fixed=TRUE);
    pattern <- "[ ]*([^ (]+)[ ]*(()|[(][^]]+[)])";
    pkgs <- gsub(pattern, "\\1", desc);
    vers <- gsub(pattern, "\\2", desc);
    vers <- gsub("[()]", "", vers);
    names(vers) <- pkgs;
    dups <- duplicated(pkgs);
    vers <- vers[!dups];
    minVersion <- vers["BatchJobs"];
    minVersion <- gsub("[>= ]*", "", minVersion);
  }
  .usePackage("BatchJobs", minVersion=minVersion, ...);
} # .useBatchJobs()


setMethodS3("getBatchJobRegistryIdByMethod", "default", function(class, method, version=NULL, ..., verbose=FALSE) {
  # Argument 'method':
  method <- Arguments$getCharacter(method);

  # Argument 'class':
  class <- Arguments$getCharacters(class);

  # Argument 'version':
  version <- Arguments$getCharacter(version);


  # Construct key from all object/arguments
  key <- list(
    class=class,
    method=method,
    version=version,
    ...
  );
  keyId <- R.cache::getChecksum(key, algo="crc32");

  id <- sprintf("%s_%s_%s", key$method, key$class[1L], keyId);

  id;
}, protected=TRUE)


setMethodS3("getBatchJobRegistryIdByMethod", "GenericDataFileSet", function(object, ...) {
  files <- getFiles(object);
  keys <- lapply(files, FUN=function(file) {
    list(filename=getFilename(file), fileSize=getFileSize(file));
  });
  getBatchJobRegistryIdByMethod(class=class(object), fileKeys=keys, ...);
}, protected=TRUE)



setMethodS3("getBatchJobRegistry", "default", function(..., skip=TRUE) {
  .useBatchJobs();
  require("BatchJobs"); # To please R CMD check.

  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'skip':
  skip <- Arguments$getLogical(skip);


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Setup
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Constructor BatchJobs registry ID...
  id <- getBatchJobRegistryIdByMethod(...);

  # ...and registry directory
  rootPath <- ".batchJobsRegistries";
  path <- file.path(rootPath, id);
  path <- Arguments$getWritablePath(path);

  # Allow comma-separated registry directory paths
  oopts <- options("BatchJobs.check.posix");
  on.exit({ options(oopts) }, add=TRUE);
  options("BatchJobs.check.posix"=FALSE);

  # Load BatchJobs registry or create if missing
  packages <- c("aroma.seq");
  reg <- BatchJobs::makeRegistry(id=id, file.dir=path, skip=skip, packages=packages);

  reg;
}, protected=TRUE) # getBatchJobRegistry()



############################################################################
# HISTORY:
# 2013-08-26
# o Added getBatchJobRegistry() and getBatchJobRegistryIdByMethod().
# o Added .usePackage() and .useBatchJobs().
# o Created.
############################################################################
