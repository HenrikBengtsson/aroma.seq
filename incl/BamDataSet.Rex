filename <- "AB042N.bowtie.sorted.bam";
path <- "bamData/AlbertsonD_2012_Bladder_VUMC/Generic/";
df <- BamDataFile(filename, path=path);
print(df);

dfList <- list(df, df, df);
ds <- BamDataSet(dfList);
print(ds);
