# $1 is the JA folder
# $2 is the sample name
cutadapt -O 3 -U 1 -m 15 --nextseq-trim=20 -j 24 --trim-n -q 20 -a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	-g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-o Trim/$2.1.fq.gz -p Trim/$2.2.fq.gz $DATA/$1/$2.1.fq.gz $DATA/$1/$2.2.fq.gz 1>Trim/$2.log

