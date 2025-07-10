#  put the name of the sample in ID
#  using command such as
#  l -d Pl* | rev | cut -f 1 -d " " | rev | cut -f 1 -d "." > ID
echo -e "ID\tA\tC\tG\tT\tN" > Daniel/Results/nt_frq.txt
for i in HX{01..04}
do
	numA=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' |tr [CGTN\\n] \\t|sed 's/	//g'|wc -c)
	numC=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' |tr [AGTN\\n] \\t|sed 's/	//g'|wc -c)
	numG=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' |tr [ACTN\\n] \\t|sed 's/	//g'|wc -c)
	numT=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' |tr [ACGN\\n] \\t|sed 's/	//g'|wc -c)
	numN=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' |tr [CGTA\\n] \\t|sed 's/	//g'|wc -c)
	echo -e $i"\t"$numA"\t"$numC"\t"$numG"\t"$numT"\t"$numN >> Daniel/Results/nt_frq.txt
done
