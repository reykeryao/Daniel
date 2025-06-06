#  put the name of the sample in ID
#  using command such as
#  l -d Pl* | rev | cut -f 1 -d " " | rev | cut -f 1 -d "." > ID
echo -e "ID\tA\tC\tG\tT\tN" > nt_frq.txt
for i in HX{01..04}
do
	numA=$(awk 'NR%2==0 {print}' Merged/$i.fa |tr [CGTN\\n] \\t|sed 's/	//g'|wc -c)
	numC=$(awk 'NR%2==0 {print}' Merged/$i.fa |tr [AGTN\\n] \\t|sed 's/	//g'|wc -c)
	numG=$(awk 'NR%2==0 {print}' Merged/$i.fa |tr [ACTN\\n] \\t|sed 's/	//g'|wc -c)
	numT=$(awk 'NR%2==0 {print}' Merged/$i.fa |tr [ACGN\\n] \\t|sed 's/	//g'|wc -c)
	numN=$(awk 'NR%2==0 {print}' Merged/$i.fa |tr [CGTA\\n] \\t|sed 's/	//g'|wc -c)
	echo -e $i"\t"$numA"\t"$numC"\t"$numG"\t"$numT"\t"$numN >> nt_frq.txt
done
