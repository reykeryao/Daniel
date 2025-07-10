kraken2 --use-names --db $REF/TGSEQ/Kraken_db/SPF/ --threads 12 --gzip-compressed --output $1.1.results Trim/$1.1.fq.gz &
kraken2 --use-names --db $REF/TGSEQ/Kraken_db/SPF/ --threads 12 --gzip-compressed --output $1.2.results Trim/$1.2.fq.gz &
wait
paste $1.1.results $1.2.results | awk '($2==$7) && ($1=="U") && ($6=="U") {print $2}' FS=\\t OFS=\\t > $1.ID
seqkit grep -f $1.ID Trim/$1.1.fq.gz -o Daniel/Filtered_by_kraken/$1.1.fq.gz &
seqkit grep -f $1.ID Trim/$1.2.fq.gz -o Daniel/Filtered_by_kraken/$1.2.fq.gz &
wait
rm $1.ID
rm $1.1.results
rm $1.2.results
