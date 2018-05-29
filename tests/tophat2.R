library("aroma.seq")

if (isCapableOf(aroma.seq, "tophat2")) {
  bin <- findTopHat2()
  print(bin)
}
