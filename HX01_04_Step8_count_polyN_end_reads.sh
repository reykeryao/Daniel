#  put the name of the sample in ID
#  using command such as
#  l -d Pl* | rev | cut -f 1 -d " " | rev | cut -f 1 -d "." > ID
echo -e "ID\tTotal\tA2\tA3\tA4\tA5\tA6\tA7\tA8\tA9\tA10\tC2\tC3\tC4\tC5\tC6\tC7\tC8\tC9\tC10\tG2\tG3\tG4\tG5\tG6\tG7\tG8\tG9\tG10\tT2\tT3\tT4\tT5\tT6\tT7\tT8\tT9\tT10" > end_polyN.txt
for i in HX{01..04}
do
	Total=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | wc -l)
	A2=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "AA$" | wc -l)
	A3=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "AAA$" | wc -l)
	A4=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "AAAA$" | wc -l)
	A5=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "AAAAA$" | wc -l)
	A6=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "AAAAAA$" | wc -l)
	A7=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "AAAAAAA$" | wc -l)
	A8=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "AAAAAAAA$" | wc -l)
	A9=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "AAAAAAAAA$" | wc -l)
	A10=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "AAAAAAAAAA$" | wc -l)
	
  C2=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "CC$" | wc -l)
  C3=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "CCC$" | wc -l)
  C4=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "CCCC$" | wc -l)
  C5=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "CCCCC$" | wc -l)
  C6=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "CCCCCC$" | wc -l)
  C7=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "CCCCCCC$" | wc -l)
  C8=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "CCCCCCCC$" | wc -l)
  C9=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "CCCCCCCCC$" | wc -l)
  C10=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "CCCCCCCCCC$" | wc -l)

  G2=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "GG$" | wc -l)
  G3=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "GGG$" | wc -l)
  G4=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "GGGG$" | wc -l)
  G5=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "GGGGG$" | wc -l)
  G6=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "GGGGGG$" | wc -l)
  G7=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "GGGGGGG$" | wc -l)
  G8=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "GGGGGGGG$" | wc -l)
  G9=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "GGGGGGGGG$" | wc -l)
  G10=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "GGGGGGGGGG$" | wc -l)

  T2=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "TT$" | wc -l)
  T3=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "TTT$" | wc -l)
  T4=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "TTTT$" | wc -l)
  T5=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "TTTTT$" | wc -l)
  T6=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "TTTTTT$" | wc -l)
  T7=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "TTTTTTT$" | wc -l)
  T8=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "TTTTTTTT$" | wc -l)
  T9=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "TTTTTTTTT$" | wc -l)
  T10=$(zcat Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "TTTTTTTTTT$" | wc -l)

	echo -e $i"\t"$Total"\t"$A2"\t"$A3"\t"$A4"\t"$A5"\t"$A6"\t"$A7"\t"$A8"\t"$A9"\t"$A10"\t"$C2"\t"$C3"\t"$C4"\t"$C5"\t"$C6"\t"$C7"\t"$C8"\t"$C9"\t"$C10"\t"$G2"\t"$G3"\t"$G4"\t"$G5"\t"$G6"\t"$G7"\t"$G8"\t"$G9"\t"$G10"\t"$T2"\t"$T3"\t"$T4"\t"$T5"\t"$T6"\t"$T7"\t"$T8"\t"$T9"\t"$T10  >> end_polyN.txt
done
