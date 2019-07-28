findBowtie2 <- function(..., command=c("bowtie2", "bowtie2-build", "bowtie2-inspect", "bowtie2-align", "bowtie2-align-l", "bowtie2-align-s")) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'command':
  command <- match.arg(command)

  # Aroma-specific variable
  path <- getExternalHome("BOWTIE2_HOME")

  versionPattern <- c("-version"=".*version ([0-9.]+).*")
  findExternal(command=command, path=path, versionPattern=versionPattern, ...)
} # findBowtie2()


queryBowtie2 <- function(what=c("support:fastq.gz"), ...) {
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'what':
  what <- match.arg(what)

  bin <- findBowtie2()
  ver <- attr(bin, "version")

  res <- FALSE
  if (is.null(ver)) ver <- NA

  if (what == "support:fastq.gz") {
    perl <- findPerl()
    suppressWarnings({
    resT <- system2(perl, args='-e "use POSIX; mkfifo(\'/tmp/aroma.seq-bowtie2-query\', 0700);"', stdout=TRUE, stderr=TRUE)
    })
    resT <- paste(resT, collapse="\n")
    supported <- (regexpr("POSIX::mkfifo not implemented", resT) == -1L)
    res <- supported
    if (!supported) {
      why <- sprintf("Your bowtie2 (v%s) does not support reading gzipped FASTQ files on this platform (%s)", ver, .Platform$OS.type)
      attr(res, "why") <- why
    }
  }

  attr(res, "version") <- ver

  res
} # queryBowtie2()
