.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname);
  pkg <- AromaSeq(pkgname);
  assign(pkgname, pkg, envir=ns);
} # .onLoad()

.onAttach <- function(libname, pkgname) {
  pkg <- get(pkgname, envir=getNamespace(pkgname));

  startupMessage(pkg, '\n\n',
    '-------------------------------------------------\n',
    ' During developing phase, install/update using:\n',
    '   source("http://aroma-project.org/hbLite.R")\n',
    '   hbInstall("aroma.seq", devel=TRUE)\n',
    '-------------------------------------------------\n'
  );
} # .onAttach()
