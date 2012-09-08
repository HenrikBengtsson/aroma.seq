############################################################################
# 
#
############################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reference genome FASTA file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path=annotationData/organisms/LambdaPhage
echo "Annotation path: $path"
if ! test -d $path; then
  echo "No such path: $path"
  exit 0
fi
filename=lambda_virus.fa
pathname=$path/$filename
echo "Annotation pathname: $pathname"
if ! test -f $pathname; then
  echo "No such pathname: $pathname"
  exit 0
fi

ls -l $path


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Single-end alignment
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bwa aln $pathname fastqData/LambdaVirusExample/Generic/reads_1.fq > foo.sai 2> foo.sai.stderr



############################################################################
# HISTORY:
# 2012-09-07
# o Created.
############################################################################
