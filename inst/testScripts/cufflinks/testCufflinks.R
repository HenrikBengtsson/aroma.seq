
library(R.utils)
source("cufflinks.R")
source("findCmd.R")
source("systemCufflinks.R")

## cufflinks -o cufflinksOut accepted_hits.bam


if (FALSE) {
  cufflinks(bams="~/workspace/projects/2012/20120612.CCSP/201302/tophatLevel2/test/MyOutput/accepted_hits.bam",
            cufflinksOptions=list(o="cufflinksOut"))
} else if (TRUE) {
  systemCufflinks(commandName="cufflinks",
                  args=list("~/workspace/projects/2012/20120612.CCSP/201302/tophatLevel2/test/MyOutput/accepted_hits.bam",
                    "-o"="cufflinksOut"),
                  system2ArgsList=list(stderr="StdOut"),  ## This tests the capture of stderr; do we want this to be a default?
                  .fake=FALSE, verbose=FALSE)
}





