hisat2 -p 12 -x Ref/template -U Merged/$1.fq.gz | samtools view -f16 - > $1.sam
cut -f 10 $1.sam | sort | uniq -c | sort -k1,1nr | awk '{i++;print ">ID"i"_"$1"\n"$2}' > MAFFT/$1.fa
rm $1.sam
