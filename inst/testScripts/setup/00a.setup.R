library("R.utils");
stopifnot(packageVersion("R.utils") >= "1.16.5");

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Install extra packages, iff missing
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
if (!isPackageInstalled("R.menu")) {
  source("http://aroma-project.org/hbLite.R");
  hbLite("R.menu");
}


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup link to test scripts directory
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- "testScripts/";
path <- Arguments$getReadablePath(path, mustExist=FALSE);
if (!isDirectory(path)) {
  pathT <- system.file("testScripts", package="aroma.affymetrix");
  pathT <- Arguments$getReadablePath(pathT);
  createLink(target=pathT);
}
path <- Arguments$getReadablePath(path);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create a launch.R file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
pathname <- "launch.R";
if (!isFile(pathname)) {
  code <- c(
    '# Usage: source("launch.R")',
    'if (interactive()) savehistory();',
    'library("R.menu");',
    'launchMenu("testScripts");'
  );
  cat(file=pathname, code, sep="\n");
}
pathname <- Arguments$getReadablePathname(pathname);
