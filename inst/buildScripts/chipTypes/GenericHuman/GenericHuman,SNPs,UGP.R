if (interactive()) savehistory();
library("aroma.seq");
library("R.menu");
verbose <- Verbose(threshold=-10, timestamp=TRUE);
options(width=60);

chipType <- "GenericHuman";


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# User settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
setOption(aromaSettings, "user/initials", "HB");
setOption(aromaSettings, "user/fullname", "Henrik Bengtsson");
obf <- sprintf("%s@%s", "henrik.bengtsson", "aroma-project.org");
setOption(aromaSettings, "user/email", obf);
saveAnywhere(aromaSettings);

fullname <- getOption(aromaSettings, "user/fullname");
stopifnot(!is.null(fullname));
email <- getOption(aromaSettings, "user/email");
stopifnot(!is.null(email));
user <- getOption(aromaSettings, "user/initials");
stopifnot(!is.null(user));


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Settings
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
snpVersions <- c("snp135");
choice <- textMenu(snpVersions, title="Choose SNP version: ", value=FALSE);
datestamp <- format(Sys.Date(), format="%Y%m%d");



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup the BED file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path <- file.path("annotationData", "genomes", "Human");
snpVersionTag <- "snp135";
filename <- sprintf("%s-UCSCGenomeBrowser.bed", snpVersionTag);
db <- TabularTextFile(filename, path=path, skip=1L);
setColumnNameTranslator(db, function(...) {
  c("chromosome", "start", "stop", "unitName", "x", "strand");
});
print(db);



# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create UGP from BED file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
data <- readDataFrame(db, colClassPatterns=c("chromosome"="character", "start"="integer"));
str(data);

chr <- data$chromosome;
chr <- gsub("^chr", "", chr);
chr[chr == "X"] <- 23;
chr[chr == "Y"] <- 24;
chr[chr == "M"] <- 25;
chr <- as.integer(chr);
data$chromosome <- chr;
keep <- is.finite(chr);
data <- data[keep,];
str(data);

nbrOfUnits <- nrow(data);
chipType <- "GenericHuman";
chrsTag <- sprintf("chr%s", seqToHumanReadable(data$chromosome));
tags <- sprintf("SNPs,%s,%s,%s%s", chrsTag, snpVersionTag, user, datestamp);
ugp <- NULL;
tryCatch({
  ugp <- AromaUgpFile$byChipType(chipType, tags=tags, nbrOfUnits=nbrOfUnits);
}, error = function(ex) {})
if (is.null(ugp)) {
  path <- file.path("annotationData", "chipTypes", chipType);
  fullname <- paste(c(chipType, tags), collapse=",");
  filename <- sprintf("%s.ugp", fullname);
  ugp <- AromaUgpFile$allocate(filename, path=path, chipType=chipType, platform="Generic", nbrOfRows=nbrOfUnits);
}
print(ugp);

# Write
ugp[,1] <- data$chromosome;
ugp[,2] <- data$start;


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Update the file footer
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
ftr <- readFooter(ugp);
ftr$createdBy <- list(
  fullname = fullname, 
  email = email
);
ftr$srcFile <- list(
  filename=getFilename(db), 
  filesize=getFileSize(db), 
  checksum=getChecksum(db)
);
writeFooter(ugp, ftr);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# Statistics
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
print(ugp);

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WHAT'S NEW:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
print(table(ugp[,1], exclude=NULL));

##      1      2      3      4      5      6      7      8
## 862464 920182 778682 797371 703636 748582 636497 615317
##      9     10     11     12     13     14     15     16
## 479177 556367 543581 526795 402736 361481 324059 349745
##     17     18     19     20     21     22     23     24
## 299864 316677 247408 250148 162748 153284 286306   6629
##     25   <NA>
##    154      0
