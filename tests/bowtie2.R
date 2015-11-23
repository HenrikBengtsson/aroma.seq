library("aroma.seq")

if (isCapableOf(aroma.seq, "bowtie2")) {
  bin <- findBowtie2()
  print(bin)
  ver <- attr(bin, "version")

  bin <- findBowtie2(command="bowtie2")
  print(bin)

  bin <- findBowtie2(command="bowtie2-build")
  print(bin)

  bin <- findBowtie2(command="bowtie2-inspect")
  print(bin)

  if (!is.na(ver)) {
    if (ver >= "2.2.0") {
      ## Since Bowtie2 v2.2.0 (2014-02-17)
      bin <- findBowtie2(command="bowtie2-align-l")
    } else {
      bin <- findBowtie2(command="bowtie2-align")
    }
    print(bin)
  }
}
