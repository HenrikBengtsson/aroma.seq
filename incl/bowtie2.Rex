# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Low-level wrapper to bowtie2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

# Command line usage:
#  $BT2_HOME/bowtie2 -x lambda_virus -U $BT2_HOME/example/reads/reads_1.fq 
#                    -S eg1.sam --seed 42
# R usage alt #1:
bowtie2("-x"="lambda_virus", "-U"="$BT2_HOME/example/reads/reads_1.fq", "-S"="eg1.sam", "--seed"=42);

# R usage alt #2:
bowtie2(x="lambda_virus", U="$BT2_HOME/example/reads/reads_1.fq", S="eg1.sam", seed=42);
# Here bowtie2() should internally name the command line options as
# -x, -U, -S and --seed, by the rule that arguments with one character names
# should have a single dash prefix and the ones with two or more characters
# should have a double dash prefix.
