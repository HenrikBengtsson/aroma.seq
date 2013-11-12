## Look for existing generic functions also in imported namespaces.
## This will affect whether setGenericS3() creates a generic function
## or not.
options("R.methodsS3:checkImports:setGenericS3"=TRUE)


# setGenericS3() should always call the arguments for replacement
# functions (..., value).  The last one must be called 'value'!
`directoryStructure<-` <- function(..., value) {
  UseMethod("directoryStructure<-")
}

