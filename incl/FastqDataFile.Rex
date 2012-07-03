filename <- "reads_1.fq";
path <- "fastqData/LambdaVirusExample/Generic/";
df <- FastqDataFile(filename, path=path);
print(df);

dfList <- list(df, df, df);
ds <- FastqDataSet(dfList);
print(ds);
