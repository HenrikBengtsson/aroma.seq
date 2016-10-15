.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname)
  pkg <- AromaSeq(pkgname)
  assign(pkgname, pkg, envir=ns)
} # .onLoad()

.onAttach <- function(libname, pkgname) {
  pkg <- get(pkgname, envir=getNamespace(pkgname))

  msg <- c(
    'During developing phase, install/update using:',
    'source("http://callr.org/install#HenrikBengtsson/aroma.seq")'
  )

  # Enable automate parallel processing via futures
  setOption("R.filesets/parallel", "future")
  msg <- c(msg,
    '',
    'Parallel processing enabled using futures',
    'To disable: setOption("R.filesets/parallel", "none")'
  )

  startupMessage(pkg, '\n\n',
    '------------------------- aroma.seq -------------------------\n',
    paste(c(msg, ''), collapse="\n"),
    '-------------------------------------------------------------\n'
  )
} # .onAttach()
