pi3 <- function(x) {
  s <- prettyNum(round(x), big.mark=",")
  names(s) <- names(x)
  s
} # pi3()


