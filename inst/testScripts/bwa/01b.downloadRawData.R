path <- system.file("testScripts/R", package="aroma.affymetrix");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

library("aroma.seq");
verbose && enter(verbose, "Downloading raw data");


##########################################################################
# Data set:
# LambdaVirusExample/
#   Generic/
#    reads_1.fq reads_1.fq [2]
#
# The example data that comes with the bowtie2 software.
#
# URL: http://bowtie-bio.sourceforge.net/bowtie2/
##########################################################################
rootPath <- "fastqData";
dataSet <- "LambdaVirusExample";
platform <- "Generic";

verbose && cat(verbose, "Data set: ", dataSet);

path <- filePath(rootPath, dataSet, platform, expandLinks="any");
ds <- tryCatch({
  FastqDataSet$byPath(path);
}, error = function(ex) FastqDataSet());
if (length(ds) < 2) {
  downloadBowtie2ExampleData();
}

ds <- FastqDataSet$byPath(path);
print(ds);

verbose && exit(verbose);
