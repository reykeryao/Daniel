mafft --auto --addfragments MAFFT/$1.fa --reorder --thread -12 --keeplength --preservecase Ref/mafft.fa > MAFFT/$1.output
fasta_formatter -t -i MAFFT/$1.output -o $1.tmp
mv $1.tmp MAFFT/$1.output
