library("aroma.seq")

cat("Tools supported by the package:\n")
print(capabilitiesOf(aroma.seq))

# Check whether BWA is supported
print(isCapableOf(aroma.seq, "bwa"))
