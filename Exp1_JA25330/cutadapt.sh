# $1 is the JA folder
# $2 is the sample name
cutadapt -O 3 -m 40 -j 24 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	-A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-o Trim/$2.1.fq.gz -p Trim/$2.2.fq.gz $DATA/$1/$2.1.fq.gz $DATA/$1/$2.2.fq.gz 1>Trim/$2.log

