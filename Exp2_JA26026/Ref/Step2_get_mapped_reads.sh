hisat2 -p 24 -k 1 --norc -x Ref/template -U Merged/$1.processed.fq.gz 2> BAM/$1.log | samtools view -bS - > BAM/$1.bam
samtools view -F4 BAM/$1.bam | cut -f 10 | sort | uniq -c | sort -k1,1nr | awk '{i++;print ">ID"i"_"$1"\n"$2}' > MAFFT/$1.fa

