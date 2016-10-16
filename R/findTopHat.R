findTopHat <- function(..., command="tophat", path=NULL) {
  if (is.null(path)) {
    # Aroma-specific variable
    path <- getExternalHome("TOPHAT_HOME");
  }
  versionPattern <- c("--version"=".*TopHat[ ]*v([0-9.-]+).*");
  findExternal(command=command, path=path, versionPattern=versionPattern, ...);
} # findTopHat()


findTopHat1 <- function(..., command="tophat", version=c(1,2)) {
  # Aroma-specific variable
  path <- getExternalHome("TOPHAT1_HOME");
  findTopHat(..., command=command, path=path, version=version);
} # findTopHat1()


findTopHat2 <- function(..., command="tophat2", version=c(2,3)) {
  # Aroma-specific variable
  path <- getExternalHome("TOPHAT2_HOME");
  res <- tryCatch({
    findTopHat(..., command=command, path=path, version=version);
  }, error = function(ex) { NULL });
  if (is.null(res)) {
    res <- findTopHat(..., command="tophat", version=version);
  }
  res;
} # findTopHat2()


############################################################################
# HISTORY:
# 2013-04-01
# o CLEANUP: Now findTopHat1() and  findTopHat2() utilizes findTopHat().
# o Now findTopHat() cached the results.
# o Renamed from findTopHatv() to findTopHat().
# 2013-01-24
# o Created.
############################################################################
