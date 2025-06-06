for i in HX{01..10};
	do
		cutadapt -O 3 -U 1 -m 15 --nextseq-trim=20 -j 24 --trim-n -q 20 -a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
			-g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
			-o Trim/$i.1.fq.gz -p Trim/$i.2.fq.gz $DATA/JA25159/$i.1.fq.gz $DATA/JA25159/$i.2.fq.gz 1>Trim/$i.log
		done

