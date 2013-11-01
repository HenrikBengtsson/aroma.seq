# Run jobs in parallel on local machine
cluster.functions <- makeClusterFunctionsMulticore(
  ncpus=parallel::detectCores()
)

