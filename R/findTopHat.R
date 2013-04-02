findTopHat <- function(..., command="tophat") {
  versionPattern <- c("--version"=".*TopHat[ ]*v([0-9.-]+).*");
  findExternal(command=command, versionPattern=versionPattern, ...);
} # findTopHat()


findTopHat1 <- function(...) {
  findTopHat(..., command="tophat", version=1);
} # findTopHat1()


findTopHat2 <- function(...) {
  res <- tryCatch({
    findTopHat(..., command="tophat", version="2");
  }, error = function(ex) { NULL });
  if (is.null(res)) {
    res <- findTopHat(..., command="tophat2", version=2);
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
