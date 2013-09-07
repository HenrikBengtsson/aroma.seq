# proto.R
# Pipeline from http://www.cureffi.org/2012/09/07/an-alternative-exome-sequencing-pipeline-using-bowtie2-and-samtools/


## Look up Bioconductor varianttools



# Barebones outline, from https://sites.google.com/site/ccsp673/home/dna-seq/exome-seq/example-pipeline:
# fastqc
# bowtie2
# samtools
# annovar
# vcftools
# bedtools coverageBed
# Python PostgreSQL

# Step 0.1
# wget ftp://ftp.cbcb.umd.edu/pub/data/bowtie2_indexes/incl/hg19.zip
# unzip hg19.zip

# Step 0.2
# wget ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
# gunzip chromFa.tar.gz
# cat chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chr20.fa chr21.fa chr22.fa chrX.fa chrY.fa chrM.fa chr1_gl000191_random.fa chr1_gl000192_random.fa chr4_gl000193_random.fa chr4_gl000194_random.fa chr7_gl000195_random.fa chr8_gl000196_random.fa chr8_gl000197_random.fa chr9_gl000198_random.fa chr9_gl000199_random.fa chr9_gl000200_random.fa chr9_gl000201_random.fa chr11_gl000202_random.fa chr17_gl000203_random.fa chr17_gl000204_random.fa chr17_gl000205_random.fa chr17_gl000206_random.fa chr18_gl000207_random.fa chr19_gl000208_random.fa chr19_gl000209_random.fa chr21_gl000210_random.fa chrUn_gl000211.fa chrUn_gl000212.fa chrUn_gl000213.fa chrUn_gl000214.fa chrUn_gl000215.fa chrUn_gl000216.fa chrUn_gl000217.fa chrUn_gl000218.fa chrUn_gl000219.fa chrUn_gl000220.fa chrUn_gl000221.fa chrUn_gl000222.fa chrUn_gl000223.fa chrUn_gl000224.fa chrUn_gl000225.fa chrUn_gl000226.fa chrUn_gl000227.fa chrUn_gl000228.fa chrUn_gl000229.fa chrUn_gl000230.fa chrUn_gl000231.fa chrUn_gl000232.fa chrUn_gl000233.fa chrUn_gl000234.fa chrUn_gl000235.fa chrUn_gl000236.fa chrUn_gl000237.fa chrUn_gl000238.fa chrUn_gl000239.fa chrUn_gl000240.fa chrUn_gl000241.fa chrUn_gl000242.fa chrUn_gl000243.fa chrUn_gl000244.fa chrUn_gl000245.fa chrUn_gl000246.fa chrUn_gl000247.fa chrUn_gl000248.fa chrUn_gl000249.fa > hg19-bt2.fa
# samtools faidx hg19-bt2.fa

# Step 1 - decompress .fq files
# gunzip -c 1_1.fq.gz 1_1.fq

# Step 2 - QC
# fastqc 1_1.fq -f fastq -o fastqc/

# Step 3 - alignment
# bowtie2 --end-to-end --very-fast --rg-id "@RG\tID:1\tLB:project_name\tSM:1\tPL:ILLUMINA" -x hg19 -q -1 1_1.fq -2 2_2.fq | samtools view - -Sb -h -t hg19-bt2.fa.fai -o 1.bam

# Step 4: sort
# samtools sort 1.bam 1.srtd -m 8000000000

# Step 5: reheader (!)
# samtools view -H 1.srtd.bam > 1originalheader.sam
# samtools reheader 1reheader.sam 1.srtd.bam > 1.srtd.reh.bam

# Step 6: remove PCR duplicates
# samtools rmdup 1.srtd.reh.bam 1.srtd.reh.ddup.bam

# Step 7: index
# samtools index 1.srtd.reh.ddup.bam

# Step 8: call variants
# ls *.srtd.reh.ddup.bam > bamlist.txt
# samtools mpileup -d 200 -D -B -f ../hg19/fasta/hg19-bt2.fa -b bamlist.txt -l ../bed/SeqCap_EZ_Exome_v2.bed -u | bcftools view - -v -c -g > variants.vcf

# Step 9: annotating variants: annovar
# - Download stuff, and eventually 'convert your vcf to an annovar file'

# Step 10: converting to plink format: vcftools

# Step 11: calculate coverage: bedtools
# - Use coverageBed...

# Step 12: converting results to database format: Python and PostgreSQL











  










