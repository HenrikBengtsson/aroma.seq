############################################################################
# 
#
# REFERENCES:
# [1] Getting started with Bowtie 2: Lambda phage example, 2012-08-31,
#     http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
############################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Session information
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
echo "BT2_HOME: $BT2_HOME"


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Indexing a reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$BT2_HOME/bowtie2-build $BT2_HOME/example/reference/lambda_virus.fa lambda_virus


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Paired-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$BT2_HOME/bowtie2 -x lambda_virus -1 $BT2_HOME/example/reads/reads_1.fq -2 $BT2_HOME/example/reads/reads_2.fq -S eg2.sam

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Results
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
head eg1.sam


############################################################################
# HISTORY:
# 2012-08-31
# o Created.
############################################################################
