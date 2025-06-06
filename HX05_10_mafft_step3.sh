cd MAFFT
for i in HX{05..10}
do
        mafft --auto --addfragments $i.fa --reorder --thread -1 --keeplength --preservecase ../Ref/mafft.fa > $i".output"
        fasta_formatter -t -i $i".output" -o $i.tmp
        mv $i.tmp $i".output"
	
done &

wait
echo "all done!"
cd ..
