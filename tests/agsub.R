library("aroma.seq")

dataset <- "MyDataSet"
organism <- "HomoSapiens"
samples <- c(A="LH001", B="LH003")
samplesT <- c(sprintf("Sample_%s", samples), samples)
x <- file.path(dataset, organism, samplesT)
x[3] <- gsub("/", "_", x[3], fixed=TRUE)
names(x) <- names(samples)
print(x)

pattern <- "(.*)/(.*)/(Sample_|)([^,]*)"

replacements <- list(
  A=list(dataset="\\1", organism="\\2", sample="\\4"),
  B=c(dataset="\\1", organism="\\2", sample="\\4"),
  C=c("\\1", "\\2", "\\4")
)

for (as in c("list", "data.frame", "matrix")) {
  for (rr in seq_along(replacements)) {
    replacement <- replacements[[rr]];
    message(sprintf("\nas='%s', replacement='%s':", as, names(replacements)[rr]))
    res <- agsub(pattern, replacement, x, as=as);
    print(res)

    # Validate
    if (as == "list") {
      stopifnot(all(is.na(res[[1L]]) | res[[1L]] == dataset))
      stopifnot(all(is.na(res[[2L]]) | res[[2L]] == organism))
      stopifnot(all(is.na(res[[3L]]) | res[[3L]] == samples))
    } else {
      stopifnot(all(is.na(res[,1L]) | res[,1L] == dataset))
      stopifnot(all(is.na(res[,2L]) | res[,2L] == organism))
      stopifnot(all(is.na(res[,3L]) | res[,3L] == samples))
    }
  } # for (rr ...)
} # for (as ...)
