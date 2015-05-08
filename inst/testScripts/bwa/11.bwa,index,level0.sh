############################################################################
# 
#
############################################################################

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Reference genome FASTA file
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
path=annotationData/organisms/Lambda_phage
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

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Create an BWA index of the reference genome
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bwa index -a is $pathname

ls -l $path

############################################################################
# HISTORY:
# 2012-09-07
# o Created.
############################################################################
