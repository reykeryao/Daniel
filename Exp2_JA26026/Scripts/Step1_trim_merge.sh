fastp -w 12 -g -m -c --overlap_len_require 10 -F 1 -l 15 \
--adapter_sequence=AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
--adapter_sequence_r2=AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
--merged_out Merged/$1.fq.gz -i $DATA/$1/$2.1.fq.gz -I $DATA/$1/$2.2.fq.gz \
-o Merged/$2.un.1.fq.gz -O Merged/$2.un.2.fq.gz -j Merged/$2.json -h Merged/$2.html 
