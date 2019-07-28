# Config checks, assertions, etc.

checkConfig <- function(config)
{
  
  ## Download reference fasta if required
  if (length(findFiles(path=pathLocalAnnots, pattern="[.](fa|fasta)$")) < 1)
  {
    # Get list of reference URLs from the package
    source(file.path(pathExData, "ReferenceGenomes.R"))
    RefUrls <- RefUrlsList[[config$organism]]
    sapply(RefUrls, function(loc)
    {
      fnameGZ <- basename(loc)
      # [ Following is cludgey on two counts: assuming .gz suffix; assuming gzip'd in the first place ]
      fname <- sub(".gz", "", fnameGZ)
      # TODO:  Make this agnostic to .gz or .bz2 or uncompressed
      
      # Download file if it has not been downloaded already
      if (!isFile(file.path(pathLocalAnnots, fnameGZ)) && !isFile(file.path(pathLocalAnnots, fname)))
      {
        pathLocalAnnots <- Arguments$getWritablePath(pathLocalAnnots)
        downloadFile(loc, path=".")
        # gunzip(fnameGZ)  ## most indexers can handle .gz input, so skip this
        renameFile(fnameGZ, file.path(pathLocalAnnots, fnameGZ))
      }
    })
  }
  
  # Some assertions here, e.g.
  # stopifnot("ref index exists")
  # stopifnot("gene model file exists")
  
  # bowtie2 is required by TopHat
  stopifnot(isCapableOf(aroma.seq, "bowtie2"))
  
}

