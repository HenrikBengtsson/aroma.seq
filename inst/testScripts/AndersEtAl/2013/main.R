## source("main.R")
## - Workflow from Anders et al, http://arxiv.org/pdf/1302.3685v2.pdf

## Known issues / questions:
## 1. Currently changing dirs all over the place.  Sometimes necessary since output files written there.
##    Alternatives?  Create files locally then move?  Create files locally then symbolically link??
## 2. Do we want e.g. TopHat, Sam, edgeR output to all go in the analysis dir?
##    [TAT:  I'd say yes.  TopHat is time consuming,and could go with the data, but then we may want to rerun
##     TopHat w/ different parameters, so we shouldn't assume there will be only one TopHat run associated with the data]
##    CURRENTLY THE TOPHAT OUTPUT GOES UNDER THE DATA DIR THOUGH - NEED TO THINK ABOUT THIS FURTHER.

############################################
## Set global vars

## Dir with input read fastq files
DirData <- "/compbio/data/SRA/GSE18508/"

## Root dir for reference genome indices
DirRefIndexRoot <- DirData

## Root dir for analysis
DirAnalysisRoot <- DirData

## Dir for reference genome sequence fasta file
DirRefFa <- DirData

## Reference sequence
RefFaFiles <- c("Drosophila_melanogaster.BDGP5.70.dna.toplevel.fa")

## The following two vars specify the path to the bowtie2 reference index
## [... CURRENTLY THE REF INDEX WILL GO INTO THE REF FA DIR ...]
Organism <-"DM"
IndexFilePrefix <- "BDGP5"

## Do QA of input reads?
bQAFastq <- FALSE

## Transcript (e.g.) model
## [... CURRENTLY THIS SHOULD BE IN THE SAME DIR AS THE REFERENCE FA ...]
GeneModelFile <- "Drosophila_melanogaster.BDGP5.70.gtf"

## NB:  'samples.RData' below is representative of a general metadata structure that needs to be manually created

##  [... Get rolling ... ]

############################################
## Mishmash of pre-processing steps:
## Download R/BioC packages, build ref index, fastq QA, create metadata

source("preProcess.R")

############################################
## Align the reads to reference genome using TopHat

source("doTopHat.R")


######################################################################
## Organize, sort and index the BAM files and create SAM files

source("doSam.R")


######################################################################
## Count reads

source("doCounts.R")
## [... CURRENTLY THIS USES HTSEQ-COUNT, BUT ALTERNATIVES (E.G. SUMMARIZEOVERLAPS) SHOULD BE AVAILABLE ...]


######################################################################
## Differential expression analysis with edgeR

source("doEdgeR.R")

