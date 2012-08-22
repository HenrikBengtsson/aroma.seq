## source("test.R")
# - This is for low-level testing only

## Tip:  If using emacs, start it up using
## 'bash -c emacs &' from a cmd line
## Usage:
## library(aroma.seq)
## testSystemBowtie2Build()  should work

testSystemBowtie2Build <- function(
                          OutPath = "./tmp",
                          ...
                          )
{
    library(R.oo)
    library(R.filesets)
    library(R.utils)
    ## source("systemBowtie2Build.R")

    ## Hack to get
    Bt2Bin <- Sys.which("bowtie2")
    Bt2Home <- getParent(Bt2Bin)
    InPathName <- filePath(Bt2Home, "example/reference/lambda_virus.fa", sep="")

    ## In the real aroma workflow, the following should already be done in bowtie2Build.R
    ##  (i.e. "outPath <- Arguments$getWritablePath(outPath);")
    if (!file.exists(OutPath))
    {
        mkdirs(OutPath)
    }
    IndexPrefix <- "lambda_virus"

    systemBowtie2Build(inPathname=InPathName, outPath=OutPath, name=IndexPrefix)
}



