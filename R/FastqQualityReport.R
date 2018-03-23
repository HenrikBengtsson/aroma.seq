setConstructorS3("FastqQualityReport", function(dataSet=NULL, ...) {
  # Validate arguments
  if (!is.null(dataSet)) {
    # Argument 'dataSet':
    dataSet <- Arguments$getInstanceOf(dataSet, "FastqDataSet");
  } # if (!is.null(dataSet))

  # Arguments '...':
  args <- list(...);


  extend(Object(), c("FastqQualityReport"),
    .ds = dataSet
  );
})
