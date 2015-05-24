setConstructorS3("CnvKitCopyNumberRatioFile", function(..., .verify=FALSE) {
  extend(TabularTextFile(..., .verify=.verify), "CnvKitCopyNumberRatioFile")
})

setMethodS3("getDefaultColumnClassPatterns", "CnvKitCopyNumberRatioFile", function(this, ...) {
  c("*"="NULL", chromosome="character", "(start|end)"="integer", "(log2|weight)"="numeric")
}, protected=TRUE)


## Doesn't work, because 'flavor' is passed along all the way to read.table(),
## which gives an error on 'unused argument (flavor = "PSCBS")'.   See also
## R.methods Issue #6 (https://github.com/HenrikBengtsson/R.methodsS3/issues/6)
## setMethodS3("readDataFrame", "CnvKitCopyNumberRatioFile", function(this, ..., flavor=c("asis", "PSCBS")) {
##   # Argument 'flavor':
##   flavor <- match.arg(flavor)
##
##   data <- NextMethod("readDataFrame", object=this)
##   data
## })

setMethodS3("readDataFrameForPSCBS", "CnvKitCopyNumberRatioFile", function(this, ..., ploidy=2) {
  data <- readDataFrame(this, ...)
  chr <- as.integer(data$chromosome)
  chr[is.element(data$chromosome, "X")] <- 23L
  chr[is.element(data$chromosome, "Y")] <- 24L
  chr[is.element(data$chromosome, c("M", "MT"))] <- 25L
  data.frame(chromosome=chr, x=(data$start + data$end)/2, y=ploidy*2^data$log2, w=data$weight)
})


setMethodS3("segmentByCBS", "CnvKitCopyNumberRatioFile", function(this, ..., minLength=2e6, rows=NULL, verbose=FALSE) {
  use("PSCBS")
  gapsToSegments <- PSCBS::gapsToSegments
  `sampleName<-` <- PSCBS::`sampleName<-`


  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Validate arguments
  # - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  # Argument 'minLength':
  minLength <- Arguments$getNumeric(minLength, range=c(0,Inf))

  # Argument 'verbose':
  verbose <- Arguments$getVerbose(verbose);
  if (verbose) {
    pushState(verbose);
    on.exit(popState(verbose));
  }


  verbose && enterf(verbose, "Paired PSCBS on %s", class(this)[1])
  sample <- getName(this)
  verbose && cat(verbose, "Sample name: ", sample)

  verbose && enter(verbose, "Reading data")
  data <- readDataFrameForPSCBS(this, rows=NULL)
  verbose && str(verbose, data)
  verbose && exit(verbose)

  if (minLength > 0) {
    gaps <- findLargeGaps(data, minLength=minLength)
    knownSegments <- gapsToSegments(gaps)
  } else {
    knownSegments <- NULL
  }

  fit <- segmentByCBS(data, knownSegments=knownSegments, ..., verbose=less(verbose))
  sampleName(fit) <- sample
  verbose && print(verbose, fit)
  verbose && exit(verbose)

  fit
})


setConstructorS3("CnvKitCopyNumberRatioFileSet", function(...) {
  extend(TabularTextFileSet(...), "CnvKitCopyNumberRatioFileSet")
})

setMethodS3("findByName", "CnvKitCopyNumberRatioFileSet", function(static, ..., pattern="[.]cnr$", paths=c("cnvkitData")) {
  NextMethod("findByName", paths=paths, pattern=pattern)
}, static=TRUE)

setMethodS3("byName", "CnvKitCopyNumberRatioFileSet", function(static, ..., paths=c("cnvkitData")) {
  NextMethod("byName", paths=paths)
}, static=TRUE)
