mafft --auto --addfragments Daniel/MAFFT/$1.fa --reorder --thread -12 --keeplength --preservecase Daniel/Ref/mafft.fa > Daniel/MAFFT/$1.output
fasta_formatter -t -i Daniel/MAFFT/$1.output -o $1.tmp
mv $1.tmp Daniel/MAFFT/$1.output