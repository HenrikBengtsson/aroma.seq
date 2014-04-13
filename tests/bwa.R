library("aroma.seq")

if (isCapableOf(aroma.seq, "bwa")) {
  bin <- findBWA()
  print(bin)
}
