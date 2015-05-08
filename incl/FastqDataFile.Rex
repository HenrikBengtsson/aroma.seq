filename <- "reads_1.fq"
path <- "fastqData/LambdaVirusExample/Lambda_phage/"
df <- FastqDataFile(filename, path=path)
print(df)

dfList <- list(df, df, df)
ds <- FastqDataSet(dfList)
print(ds)
