.onLoad <- function(libname, pkgname) {
  ns <- getNamespace(pkgname);
  pkg <- AromaSeq(pkgname);
  assign(pkgname, pkg, envir=ns);
} # .onLoad()

.onAttach <- function(libname, pkgname) {
  pkg <- get(pkgname, envir=getNamespace(pkgname));

  msg <- c(
    'During developing phase, install/update using:',
    'source("http://aroma-project.org/R/install#aroma.seq")'
  );

  # Enable automate parallel processing via BatchJobs?
  if (queryRCmdCheck() == "notRunning") {
    if (utils::file_test("-f", ".BatchJobs.R")) {
      setOption("R.filesets/parallel", "BiocParallel::BatchJobs");
      msg <- c(msg,
        '',
        'Parallel processing enabled (via \'./.BatchJobs.R\')',
        'To disable: setOption("R.filesets/parallel", "none")'
      );
    }
  }

  startupMessage(pkg, '\n\n',
    '---------------------- aroma.seq ----------------------\n',
    paste(c(msg, ''), collapse="\n"),
    '-------------------------------------------------------\n'
  );
} # .onAttach()
