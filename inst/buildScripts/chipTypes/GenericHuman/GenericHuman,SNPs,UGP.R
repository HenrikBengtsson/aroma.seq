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

# Drop non-interesting "chromosomes"
keep <- is.finite(chr);
data <- data[keep,];

# Drop duplicates
dups <- duplicated(data);
if (any(dups)) {
  data <- data[!dups,];
}

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
## AromaUgpFile:
## Name: GenericHuman
## Tags: SNPs,chr1-25,snp135,HB20121031
## Full name: GenericHuman,SNPs,chr1-25,snp135,HB20121031
## Pathname: annotationData/chipTypes/GenericHuman/GenericHuman,SNPs,chr1-25,snp135
## ,HB20121031.ugp
## File size: 53.87 MB (56483592 bytes)
## RAM: 0.00 MB
## Number of data rows: 11296623
## File format: v1
## Dimensions: 11296623x2
## Column classes: integer, integer
## Number of bytes per column: 1, 4
## Footer: <createdOn>20121031 14:21:24 PDT</createdOn><platform>Generic</platform>
## <chipType>GenericHuman</chipType><createdBy><fullname>GenericHuman,SNPs,chr1-25,
## snp135,HB20121031</fullname><email>henrik.bengtsson@aroma-project.org</email></c
## reatedBy><srcFile><filename>snp135-UCSCGenomeBrowser.bed</filename><filesize>445
## 617891</filesize><checksum>0e52557c814e0e1fb93482e9b7894c1e</checksum></srcFile>
## Chip type: GenericHuman
## Platform: Generic

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
# WHAT'S NEW:
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
print(table(ugp[,1], exclude=NULL));

##      1      2      3      4      5      6      7      8
## 861880 919552 778274 795895 703201 724540 636155 614952
##      9     10     11     12     13     14     15     16
## 478897 556050 543244 526533 402432 361286 323849 349564
##     17     18     19     20     21     22     23     24
## 297842 316472 247287 250035 162652 153139 286110   6628
##     25   <NA>
##    154      0
