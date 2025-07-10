hisat2 -p 12 -x Daniel/Ref/Unknow_seq -U Daniel/Merged/$1.fq.gz | samtools view -f4 - | cut -f 10 | sort | awk '{i++;print ">ID"i"\n"$1}' > Daniel/Merged/$1.fa
fasta_formatter -t -i Daniel/Merged/$1.fa -o Daniel/Merged/$1.tab
gzip Daniel/Merged/$1.fa
gzip Daniel/Merged/$1.tab
