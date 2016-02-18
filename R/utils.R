htable <- function(x, useNA="ifany", ..., fmtstr="%s [%d]") {
  t <- table(x, useNA=useNA)
  hpaste(sprintf(fmtstr, names(t), t), ...)
}

  
