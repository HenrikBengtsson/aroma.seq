## source("testS0.R")

if (TRUE) {  ## (for standalone testing)
    library(R.filesets)
    source("findTopHat.R")
}

source("systemTophatBuild.R")

system.time({
    systemTopHatBuild('-o'="testBuild",
                      '-G'="/home/tokuyasu/projects/aroma.seq/20121109.Quigley/genes_MM9_iGenomes_no_random.gtf",
                      '--transcriptome-index'= "testBuildIdx/known",
                      '-p'=12,
                      '/home/tokuyasu/projects/aroma.seq/20121109.Quigley/mm9',
                      'reads_1.fq',
                      'reads_2.fq',
                      verbose=TRUE)
})


############################################################################
# HISTORY:
# 2012-09-24
# o Created.
############################################################################


