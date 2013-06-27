library("aroma.seq")

if (isCapableOf(aroma.seq, "bowtie2")) {
  bin <- findBowtie2()
  print(bin)

  bin <- findBowtie2(command="bowtie2")
  print(bin)

  bin <- findBowtie2(command="bowtie2-align")
  print(bin)

  bin <- findBowtie2(command="bowtie2-build")
  print(bin)

  bin <- findBowtie2(command="bowtie2-inspect")
  print(bin)
}
