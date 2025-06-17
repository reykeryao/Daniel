hisat2 -p 12 -x Ref/Unknow_seq -U Merged/$1.fq.gz | samtools view -f4 - | cut -f 10 | sort | awk '{i++;print ">ID"i"\n"$1}' > Merged/$1.fa
fasta_formatter -t -i Merged/$1.fa -o Merged/$1.tab
gzip Merged/$1.fa
gzip Merged/$1.tab
