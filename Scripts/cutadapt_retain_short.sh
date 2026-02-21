# $1 is the JA folder
# $2 is the sample name
cutadapt -O 3 -U 1 -m 1 -M 15 --nextseq-trim=20 -j 24 --trim-n -q 20 -a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
	-g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
	-o Trim_1to15/$2.1.fq.gz -p Trim_1to15/$2.2.fq.gz $DATA/$1/$2.1.fq.gz $DATA/$1/$2.2.fq.gz 1>Trim_1to15/$2.log

