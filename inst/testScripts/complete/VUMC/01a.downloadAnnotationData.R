library("aroma.core");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);
ar <- AromaRepository(verbose=TRUE);

verbose && enter(verbose, "Downloading annotation data");

chipType <- "GenericHuman";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadUGP(ar, chipType, tags=c("50kb", ".*"));
verbose && cat(verbose, "UGP: ", pathname);

pathname <- downloadChipTypeFile(ar, chipType, tags=c("50kb", ".*"), ext="unc");
verbose && cat(verbose, "UNC: ", pathname);

pathname <- downloadChipTypeFile(ar, chipType, tags=c("SNPs", "chr1-25", ".*"), ext="ugp");
verbose && cat(verbose, "SNP UGP: ", pathname);

verbose && exit(verbose);
