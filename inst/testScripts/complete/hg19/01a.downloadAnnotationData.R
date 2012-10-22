path <- system.file("testScripts/R", package="aroma.seq");
pathname <- file.path(path, "downloadUtils.R");
source(pathname);

library("aroma.seq");
verbose <- Arguments$getVerbose(-8, timestamp=TRUE);

ar <- AromaRepository(verbose=TRUE);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Local functions
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
verbose && enter(verbose, "Downloading annotation data");

chipType <- "GenericHuman";
verbose && cat(verbose, "Chip type: ", chipType);

pathname <- downloadUGP(ar, chipType, tags=c("50kb", ".*"));
verbose && cat(verbose, "UGP: ", pathname);


fa <- download1000GenomesHumanReferenceFile();
print(fa);
## FastaReferenceFile:
## Name: human_g1k_v37
## Tags:
## Full name: human_g1k_v37
## Pathname: annotationData/organisms/Human/human_g1k_v37.fasta
## File size: 2.94 GB (3153506519 bytes)
## RAM: 0.00 MB
## Total sequence length: NA
## Number of sequences: 84
## Sequence names: [84] 1 dna:chromosome chromosome:GRCh37:1:1:249250621:1, 
## 2 dna:chromosome chromosome:GRCh37:2:1:243199373:1, 3 dna:chromosome
## chromosome:GRCh37:3:1:198022430:1, ..., GL000192.1 dna:supercontig
## supercontig::GL000192.1:1:547496:1


verbose && exit(verbose);
