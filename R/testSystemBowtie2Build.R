## source("test.R")

## Tip:  If using emacs, start it up using
## 'bash -c emacs &' from a cmd line


testSystemBowtie2Build <- function(
                          OutPath = "./tmp",
                          ...
                          )
{
    library(R.oo)
    library(R.filesets)
    library(R.utils)
    ## source("systemBowtie2Build.R")

    args <- list(...)

    Bt2Bin <- Sys.which("bowtie2")
    Bt2Home <- getParent(Bt2Bin)
    InPathName <- filePath(Bt2Home, "example/reference/lambda_virus.fa", sep="")
    if (!file.exists(OutPath))
    {
        mkdirs(OutPath)
    }
    IndexPrefix <- "lambda_virus"

    systemBowtie2Build(inPathname=InPathName, outPath=OutPath, name=IndexPrefix)
}



