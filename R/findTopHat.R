findTopHat <- function(..., command="tophat") {
  versionPattern <- c("--version"=".*TopHat[ ]*v([0-9.-]+).*");
  findExternal(command=command, versionPattern=versionPattern, ...);
} # findTopHat()


findTopHat1 <- function(..., command="tophat", version=c(1,2)) {
  findTopHat(..., command=command, version=version);
} # findTopHat1()


findTopHat2 <- function(..., command="tophat2", version=c(2,3)) {
  res <- tryCatch({
    findTopHat(..., command=command, version=version);
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
