## source("doEdgeR.R")

## Step 1 edgeR Create container for count data and filter features
## Load the edgeR package and use the utility function, readDGE, to read in the COUNT
## files created from htseq-count:
library("edgeR")
counts <- readDGE(samples$countf)$counts
## In edgeR, we recommend removing features without at least 1 read per million in n of the samples, where n is the size of the smallest group of replicates (here, n <- 3 for the Knockdown group). Filter these as well as non-informative (e. g., non-aligned) features using a command like:
noint <- rownames(counts) %in%
c("no_feature","ambiguous","too_low_aQual",
  "not_aligned","alignment_not_unique")
cpms <- cpm(counts)
keep <- rowSums(cpms>1)>=3 & !noint
counts <- counts[keep,]
## Visualize and inspect the count table using:
colnames(counts) <- samples$shortname
head( counts[,order(samples$condition)], 5 )
## Create a DGEList object (s container for RNA-seq count data), as follows:
d <- DGEList(counts=counts, group=samples$condition)
######################
## Step 2 edgeR Estimate normalization factors
## Estimate normalization factors using:
d <- calcNormFactors(d)
d$samples
######################
## Step 3 edgeR Inspect sample relations
## Use the plotMDS function to create a count-specific multidimensional scaling plot (shown
## in Figure 4A):
plotMDS(d, labels=samples$shortname,
        col=c("darkgreen","blue")[factor(samples$condition)])
######################
## Step 4 edgeR Estimate dispersion and conduct statistical tests according to A if a
## simple two-group comparison, or B if a complex design.
## A. Perform statistical calculations for a simple two-group comparison
## (i) Estimate  edgeR)
## For simple designs, estimate tagwise dispersion estimates using:
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)
## Create a visual representation of the mean-variance relationship using plotMeanVar
## (shown in Figure 6A), as follows:
plotMeanVar(d, show.tagwise.vars=TRUE, NBline=TRUE)
## and plotBCV (Figure 6B), as follows:
plotBCV(d)
## (ii) Test for differential  edgeR)
## [See also: TROUBLESHOOTING]
##For a simple two-group design, perform an exact test for the difference in expression
## between the two conditions:
de <- exactTest(d, pair=c("CTL","KD"))
######################
## B. Perform statistical calculations for a complex design
## (i) Estimate  edgeR)
## For more complex designs, create a design matrix to specify the factors that are expected
## to affect expression levels:
design <- model.matrix( ~ LibraryLayout + condition, samples)
design
## Estimate dispersion values, relative to the design matrix, using the Cox-Reid (CR) ad-
## justed likelihood 10, 65 , as follows:
d2 <- estimateGLMTrendedDisp(d, design)
d2 <- estimateGLMTagwiseDisp(d2, design)
## (ii) Test for differential  edgeR)
## Given the design matrix and dispersion estimates, fit the GLM to the data:
f <- glmFit(d2, design)
## Perform a likelihood ratio test, specifying the difference of interest (here, Knockdown
## versus Control, corresponding to the 3rd column of the design matrix):
lrt <- glmLRT(f, coef=3)
######################
## Step 5 edgeR Inspect the results in graphical and tabular format
## Use the topTags function to present a tabular summary of the differential expression
## statistics (Note: topTags operates on the output of either exactTest or glmLRT, while
## only the latter is shown here):
tt <- topTags(lrt, n=nrow(d))
head(tt$table)
## Inspect the depth-adjusted reads per million some of the top differentially expressed genes:
nc <- cpm(d, normalized.lib.sizes=TRUE)
rn <- rownames(tt$table)
head(nc[rn,order(samples$condition)],5)
## Create a graphical summary, such as an M (log-fold-change) versus A (log-average-
## expression) plot 66 , here showing the genes selected as differentially expressed (with a 5%
## false discovery rate; see Figure 5A):
deg <- rn[tt$table$FDR < .05]
plotSmear(d, de.tags=deg)
######################
## Step 6 edgeR Create persistent storage of results
## Save the result table as a CSV (comma-separated values) file (alternative formats are
## possible) as follows:
write.csv( tt$table, file="toptags_edgeR.csv" )


