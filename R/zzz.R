# Allows conflicts. For more information, see library() and
# conflicts() in [R] base.
.conflicts.OK <- TRUE;

.onAttach <- function(libname, pkgname) {
  pkg <- AromaSeq(pkgname);
  assign(pkgname, pkg, pos=getPosition(pkg));

  startupMessage(pkg, '\n\n',
    '-------------------------------------------------\n',
    ' During developing phase, install/update using:\n',
    '   source("http://aroma-project.org/hbLite.R")\n',
    '   hbInstall("aroma.seq", devel=TRUE)\n',
    '-------------------------------------------------\n'
  );
} # .onAttach()
