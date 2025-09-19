fastp -A -m -c --overlap_len_require 10 --merged_out Merged/$1.fq.gz -i ../Trim/$1.1.fq.gz -I ../Trim/$1.2.fq.gz -o Merged/$1.un.1.fq.gz -O Merged/$1.un.2.fq.gz -j Merged/$1.json -h Merged/$1.html 
