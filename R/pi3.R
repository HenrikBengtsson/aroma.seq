pi3 <- function(x) {
  s <- prettyNum(as.integer(x), big.mark=",")
  names(s) <- names(x)
  s
} # pi3()


