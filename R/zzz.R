.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname)
  pkg <- AromaSeq(pkgname)
  assign(pkgname, pkg, envir = ns)
} # .onLoad()

.onAttach <- function(libname, pkgname) {
  pkg <- get(pkgname, envir = getNamespace(pkgname))

  # Enable automate parallel processing via futures
  setOption("R.filesets/parallel", "future")

  startupMessage(pkg, '\n\n',
    '------------------------- aroma.seq -------------------------\n',
    'During developing phase, install/update using:               \n',
    'source("http://callr.org/install#HenrikBengtsson/aroma.seq") \n',
    '-------------------------------------------------------------\n'
  )
}
