# source("setupDirs.R")
# Create working directory structure and fill with data

setupDirs <- function(config,
                      pathLocalAnnots=NULL,
                      pathLocalData=NULL,
                      bGunzip=FALSE,
                      ...)
{
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Copy reference and input files to local dirs
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  
  # Set up local annotation dir
  if (is.null(pathLocalAnnots)) {
    pathLocalAnnots <- file.path("annotationData", "organisms", config$organism)
    pathLocalAnnots <- Arguments$getWritablePath(pathLocalAnnots)    
  }

  # Set up local data dir
  if (is.null(pathLocalData)) {
    pathLocalData <- file.path("fastqData", config$datasetName, config$organism)
    pathLocalData <- Arguments$getWritablePath(pathLocalData)
  }
 
  # Copy reference and input data files to local dirs
  if (file.exists(config$pathRef)) {
    copyDirectory(from=config$pathRef, to=pathLocalAnnots, overwrite=config$bOverwrite)
  }
  if (file.exists(config$pathData)) {
    patternData <- "fastq[.gz]*$"
    dataFiles <- findFiles(pattern=patternData, path=config$pathData, firstOnly=FALSE)
    file.copy(from=dataFiles, to=pathLocalData, overwrite=config$bOverwrite)
  }
  # If config$pathRef or config$pathData do not exist, then the reference annotations
  # or input fastq files must already exist in the proper local dirs.
  
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Gunzip input / reference data if necessary
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if (bGunzip){
    sapply(findFiles(path=pathLocalData, pattern=".gz$", firstOnly=FALSE),
           function (f) {gunzip(f, overwrite=config$bOverwrite)})
    sapply(findFiles(path=pathLocalAnnots, pattern=".gz$", firstOnly=FALSE),
           function (f) {gunzip(f, overwrite=config$bOverwrite)})
  }

  config$pathLocalAnnots <- pathLocalAnnots
  config$pathLocalData <- pathLocalData
  
  return(config)
}

