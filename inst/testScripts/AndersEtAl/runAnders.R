# source("runAnders.R")
# Anders et al, 2013 Nature Protocol
# - An aroma.seq version of the protocol runnable entirely in R


runAnders <- function(config=NULL,
                      verbose=FALSE,
                      ...)
{
  library(aroma.seq)  # Contains Anders et al metadata files for convenience
  
  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }
  
  # Config
  if (is.null(config)) {
    config <- list()
    config$Organism <- "DM"
    config$Project <- "Anders"
    config$bRPkgInstall <- FALSE    

    # Initial setup
    config$bDownloadSra <- FALSE  # When TRUE, download SRA files
    config$bFastqDump <- FALSE    # When TRUE, run fastq-dump on SRA files
    config$bDownloadReference <- FALSE
    config$bDownloadGeneModel <- FALSE
    
    if (config$bDownloadSra) {
      config$bFastqDump <- TRUE
    }
    config$SraDir <- 'Download'  # SRA file download location
  }

  
  # - - - - - - - - - - - - - - -
  # Install Bioconductor packages
  # - - - - - - - - - - - - - - -
  if (config$bRPkgInstall) {
    source("http://www.Bioconductor.org/biocLite.R")
    biocLite("BiocUpgrade")
    biocLite( c("ShortRead","DESeq", "edgeR") )
  }
  
  
  # - - - - - - - - - - - - - - -
  # Set up local directories
  # - - - - - - - - - - - - - - -
  # The 'remote' location for annots and fastq files (i.e. will make symbolic links to these)
  AnnotPathRemote <- config$SraDir
  
  # Local dirs (following 201402 style recommendations)
  AnnotPathLocal <- file.path("annotationData", "organisms", config$Organism)
  FastqPathLocal <- file.path("fastqData", config$Project, config$Organism)
  Arguments$getWritablePath(AnnotPathLocal)
  Arguments$getWritablePath(FastqPathLocal)
  
  # Awkward way of checking whether the SRA dir is writable
  bWritable <- tryCatch(Arguments$getWritablePath(config$SraDir),
                  error = function(e) {
                    print("SRA dir not writable; FastqPathRemote will be set to FastqPathLocal")
                    FALSE 
                  })
  # If SRA dir is writable, output fastq files to it; if not, output locally
  if (bWritable!=FALSE) {
    FastqPathRemote <- config$SraDir
  } else {
    FastqPathRemote <- FastqPathLocal
  }
  
  
  # - - - - - - - - - - - - - - -
  # Download Anders example data
  # - - - - - - - - - - - - - - -
  if (config$bDownloadSra) {
    # Check for fastq-dump
    if (!isCapableOf(aroma.seq, "fastqDump")) {
      warning("fastq-dump not found; download the SRA Toolkit for your OS from
              http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software")
    }
    SraDir <- Arguments$getWritablePath(config$SraDir)
    SriFile <- system.file("exData/annotationData/AndersEtAl/SraRunInfo.csv", package="aroma.seq")
    sri <- read.csv(SriFile, stringsAsFactors=FALSE)
    keep <- grep("CG8144|Untreated-", sri$LibraryName)
    sri <- sri[keep,]  # 14 single-end and 8 paired-end files
    SraFiles <- file.path(SraDir, basename(sri$download_path))
    verbose && cat(verbose, "Downloading SRA files")
    for(i in 1:nrow(sri)) {
      f <- SraFiles[i]
      if (!file.exists(f)) {
        verbose && cat(verbose, paste("Downloading", f))
        download.file(sri$download_path[i], f)        
      }
    }
    stopifnot( all(file.exists(SraFiles)) ) # assure FTP download was successful
  }
  
  # - - - - - - - - - - - - - - -
  # Convert SRA to FASTQ format
  # - - - - - - - - - - - - - - -
  if (config$bFastqDump) {
    # Given valid SRA dir, run fastq-dump on the .sra files it contains
    if (is.character(config$SraDir) && nchar(config$SraDir)>0) {
      if (!exists("SraFiles")) {
        SraFiles <- findFiles(path=SraDir, pattern="[.]sra$", firstOnly=FALSE)
      }
      for(f in SraFiles) {
        fastqDump(f, outPath=file.path(FastqPathRemote, sub("[.]sra$", "", basename(f))))
        # - For atomicity, each sra file gets its own output folder
      }
    } else {
      stop("SRA files not found; config$SraDir not a valid directory")
    }
  }
  
  
  # - - - - - - - - - - - - - - -
  # aroma.seq step: Make symbolic links to fastq files
  # - - - - - - - - - - - - - - -
  fastqFiles <- findFiles(FastqPathRemote, pattern="[.]fastq[.gz]*$", recursive=TRUE, firstOnly=FALSE)
  FastqFilesLocal <-
    sapply(fastqFiles, function(f) {
      localFile <- file.path(FastqPathLocal, basename(f))
      createLink(link=localFile, target=f)
      if (regexpr("[.]gz$", f, ignore.case=TRUE) != -1L) {
        localFile <- gunzip(localFile, skip=TRUE) # Replaces link with *gunzipped file*
      }
      localFile
    })
  
  
  # - - - - - - - - - - - - - - -
  # Download the reference genome
  # - - - - - - - - - - - - - - -
  if (config$bDownloadReference) {
    load(system.file("exData/annotationData/RefUrlsList.RData", package="aroma.seq"))
    # $ wget ftp://ftp.ensembl.org/pub/release-70/fasta/drosophila_melanogaster/dna/Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa.gz
    # $ gunzip Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa.gz
    refFileRemote <- file.path(AnnotPathRemote, basename(RefUrlsList$DM["reference"]))
    status <- download.file(from=RefUrlsList$DM["reference"],
                  to=refFileRemote)
    if (status!=0) {stop("Reference genome download failed")}
  }
  refFileRemote <- findFiles(path=AnnotPathRemote, pattern="[.]fasta$", firstOnly=TRUE)
  # - Presume that the first fasta file in AnnotPathRemote is the reference.
  refFileLocal <- createLink(link=file.path(AnnotPathLocal, basename(refFileRemote)),
                             target=refFileRemote, skip=TRUE)
  
  
  # - - - - - - - - - - - - - - -
  # Get gene model annotations
  # - - - - - - - - - - - - - - -
  if (config$bDownloadGeneModel) {
    download.file(from=RefUrlsList$DM["geneModel"],
                  to=file.path(AnnotPathRemote, basename(RefUrlsList$DM["geneModel"])))
    # shell$  wget ftp://ftp.ensembl.org/pub/release-70/gtf/drosophila <- melanogaster/Drosophila <- melanogaster.BDGP5.70.gtf.gz
    # shell$  gunzip Drosophila <- melanogaster.BDGP5.70.gtf.gz
    # CRITICAL: 'Make sure that the gene annotation uses the same coordinate system as the reference FASTA file.'
  }
  
  #####################################################################
  
  # - - - - - - - - - - - - - - -
  # Build the reference index
  # - - - - - - - - - - - - - - -
  if (!is.null(refFileLocal)) {
    # shell$  bowtie2-build -f Drosophila <- melanogaster. BDGP5.70.dna.toplevel.fa Dme1 <- BDGP5 <- 70
    refFa <- FastaReferenceFile(refFileLocal)  
    refIndex <- buildBowtie2IndexSet(refFa, verbose=-100)
  } else {
    ##  ... SEARCH FOR reference index IN THE STD PLACE ...
    bt2File <- findFiles(path=AnnotPathLocal, pattern="[.]bt2$")
    bt2Prefix <- sub("(.rev)*[.][0-9][.]bt2$", "")
    refIndex <- Bowtie2IndexSet$byPrefix(bt2Prefix)
  }
  
  
  # - - - - - - - - - - - - - - -
  # Assess sequence quality control with ShortRead  (timing ~2 h)
  # - - - - - - - - - - - - - - -
  # setwd('working directory to where the FASTQ files are situated')
  library("ShortRead")
  fqQC <- qa(dirPath=".", pattern=".fastq$", type="fastq")
  report(fqQC, type="html", dest="fastqQAreport")
  # Use a web browser to inspect the generated HTML file (here, stored in the ‘fastqQAreport’ directory) with the
  # quality-assessment report (see ANTICIPATED RESULTS for further details)
  
  
  # - - - - - - - - - - - - - - -
  # Collect metadata of experimental design (timing <1 h)
  # - - - - - - - - - - - - - - -
  # Create a table of metadata called ‘samples’ (see ‘Constructing metadata table’ in Experimental Design)
  
  
  # - - - - - - - - - - - - - - -
  # Collapse the initial table (sri) to one row per sample:
  # - - - - - - - - - - - - - - -
  sri$LibraryName <- gsub("S2_DRSC_","",sri$LibraryName) # trim label
  samples <- unique(sri[,c("LibraryName","LibraryLayout")])
  for(i in seq <- len(nrow(samples))) {
    rw <- (sri$LibraryName==samples$LibraryName[i])
    if(samples$LibraryLayout[i]=="PAIRED") {
      samples$fastq1[i] <- paste0(sri$Run[rw],"_1.fastq",collapse=",")
      samples$fastq2[i] <- paste0(sri$Run[rw],"_2.fastq",collapse=",") } else {
        samples$fastq1[i] <- paste0(sri$Run[rw],".fastq",collapse=",")
        samples$fastq2[i] <- ""
      }
  }
  
  
  # - - - - - - - - - - - - - - -
  # 5| Add important or descriptive columns to the metadata table (
  # - - - - - - - - - - - - - - -
  samples$condition <- "CTL"
  samples$condition[grep("RNAi",samples$LibraryName)] <- "KD"
  samples$shortname <- paste(substr(samples$condition,1,2),
                             substr(samples$LibraryLayout,1,2),
                             seq <- len(nrow(samples)), sep=".")
  
  
  # - - - - - - - - - - - - - - -
  # 6| Carefully inspect (and correct, if necessary) the metadata table.
  # - - - - - - - - - - - - - - -
  samples
  
  
  
  # - - - - - - - - - - - - - - -
  # Align the reads (using tophat2) to the reference genome (timing ~45 min per sample)
  # - - - - - - - - - - - - - - -
  
  # 7| By using R string manipulation, construct the Unix commands to call tophat2.
  gf <- "Drosophila_melanogaster.BDGP5.70.gtf"
  bowind <- "Dme1_BDGP5_70"
  cmd <- with(samples, paste("tophat2 -G", gf, "-p 5 -o", LibraryName, bowind, fastq1, fastq2))
  cmd
  
  
  # 8| Run these commands (i.e., copy and paste) in a Unix terminal.
  # '$ cmd.sh'
  
  
  # - - - - - - - - - - - - - - -
  # Organize, sort and index the BAM files and create SAM files (timing ~1 h)
  # - - - - - - - - - - - - - - -
  
  # 9| Organize the BAM files into a single directory, sort and index them and create SAM files by running the following R-generated commands:
  for(i in seq <- len(nrow(samples))) {
    lib <- samples$LibraryName[i]
    ob <- file.path(lib, "accepted_hits.bam")
    # sort by name, convert to SAM for htseq-count
    cat(paste0("samtools sort -n ",ob," ",lib,"_sn"),"\n")
    cat(paste0("samtools view -o ",lib,"_sn.sam ",lib,"_sn.bam"),"\n")
    # sort by position and index for IGV
    cat(paste0("samtools sort ",ob," ",lib,"_s"),"\n")
    cat(paste0("samtools index ",lib,"_s.bam"),"\n\n")
  }
  # CRITICAL:  Users should be conscious of the disk space that may get used in these operations. In the command above, sorted-by-name SAM and BAM files (for htseq-count), as well as a sorted-by-chromosome-position BAM file (for IGV), are created for each original accepted hits.bam file. User may wish to delete (some of) these intermediate files after the steps below.
  
  
  # - - - - - - - - - - - - - - -
  # Inspect alignments with IGV ● timing <20 min
  # - - - - - - - - - - - - - - -
  # 10| Start IGV, select the D. melanogaster (dm3) genome, and then load the BAM files (with s in the filename) as well as the GTF file.
  
  # 11| Zoom in on an expressed transcript until individual reads are shown and check whether the reads align at and across exon-exon junctions, as expected, given the annotation (Fig. 3).
  
  # 12| If any positive and negative controls are known for the system under study (e.g., known differential expression), direct the IGV browser to these regions to confirm that the relative read density is different according to expectation.
  
  
  # Count reads using htseq-count (timing ~3 h)
  # 13| Add the names of the COUNT files to the metadata table and call HTSeq from the following R-generated Unix commands:
  samples$countf <- paste(samples$LibraryName, "count", sep=".")
  gf <- "Drosophila_melanogaster.BDGP5.70.gtf"
  cmd <- paste0("htseq-count -s no -a 10 ", samples$LibraryName, "_sn.sam ", gf," > ", samples$countf)
  cmd
  # CRITICAL STEP:  The option -s signifies that the data are not from a stranded protocol (this may vary by experiment) and the -a option specifies a minimum score for the alignment quality.
  
  
  ############################################
  
  # 14| For differential expression analysis with edgeR, follow option A for simple designs and option B for complex designs; for differential expression analysis with DESeq, follow option C for simple designs and option D for complex designs.
  # (A) edgeR - simple design
  # (B) edgeR - complex design
  # (C) Deseq—simple design
  # (D) Deseq—complex design
  
  # 15| As another spot check, point the IGV genome browser (with GTF and BAM files loaded) to a handful of the top differen- tially expressed genes and confirm that the counting and differential expression statistics are appropriately represented.
  
  verbose && exit(verbose);
  
}  # runAnders()

