
.optionsList2String <- function(optionsList=list())
{
  ## Return if nothing to process
  if (length(optionsList) < 1) { return(NULL) }

  ## Split into logical (boolean) and non-logical (name-value pair) options
  tmpLogical <- unlist(optionsList[sapply(optionsList, is.logical)])
  tmpNonLogical <- unlist(optionsList[!sapply(optionsList, is.logical)])  ## NULLs are auto-removed

  ## For logical options: 1) get the names of the TRUE ones, 2) add a hyphen (or two)
  tmpLogical <- tmpLogical[tmpLogical==TRUE]
  nms <- names(tmpLogical)
  optionsLogical <-
    sapply(nms, function(nm)
           {
             str <- ifelse(nchar(nm)==1,
                           paste("-", nm, sep=""),
                           paste("--", nm, sep=""))
             return(str)
           })

  ## For "non-logical" options:
  ##  1) create two vectors: names and values; 2) paste together with appropriate hypenation
  nms <- names(tmpNonLogical)
  vals <- as.vector(tmpNonLogical)  ## (Could do some sanity checking of the vals here)
  hyphenVec <-
    sapply(nms, function(nm)
           {
             ifelse(nchar(nm)==1, "-", "--")
           })
  optionsNonLogical <- paste(paste(hyphenVec, nms, sep=""), vals, sep=" ")
  ## - This works appropriately even when tmpNonLogical is NULL

  ## Combine vectors into a string and return
  cmdOptions <- paste(c(optionsLogical, optionsNonLogical), collapse=" ")
  return(cmdOptions)
}

