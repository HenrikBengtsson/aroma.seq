pi3 <- function(x) {
  s <- prettyNum(round(x), big.mark=",", preserve.width="none")
  names(s) <- names(x)
  s
} # pi3()


