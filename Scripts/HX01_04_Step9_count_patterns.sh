## number of reads start with N5
echo -e "ID\tTotal\tAC\tAG\tAT\tCG\tCT\tGT\tAA\tCC\tGG\tTT" > Daniel//Results/polyN_patterns.txt
for i in HX{01..04}
do
	Total=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' | wc -l)
	AC=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep -v "AAAA[A|C|G|T|N]*CCCC" |grep -v "CCCC[A|C|G|T|N]*AAAA" |wc -l)
	AG=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep -v "AAAA[A|C|G|T|N]*GGGG" |grep -v "GGGG[A|C|G|T|N]*AAAA" |wc -l)
	AT=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep -v "AAAA[A|C|G|T|N]*TTTT" |grep -v "TTTT[A|C|G|T|N]*AAAA" |wc -l)
	CG=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep -v "CCCC[A|C|G|T|N]*GGGG" |grep -v "GGGG[A|C|G|T|N]*CCCC" |wc -l)
	CT=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep -v "CCCC[A|C|G|T|N]*TTTT" |grep -v "TTTT[A|C|G|T|N]*CCCC" |wc -l)
	GT=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep -v "GGGG[A|C|G|T|N]*TTTT" |grep -v "TTTT[A|C|G|T|N]*GGGG" |wc -l)
	AA=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "AAAA[A|C|G|T|N]*AAAA" |wc -l)
	CC=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "CCCC[A|C|G|T|N]*CCCC" |wc -l)
	GG=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "GGGG[A|C|G|T|N]*GGGG" |wc -l)
	TT=$(zcat Daniel/Merged/$i.fa.gz | awk 'NR%2==0 {print}' | grep "TTTT[A|C|G|T|N]*TTTT" |wc -l)
	
echo -e $i"\t"$Total"\t"$AC"\t"$AG"\t"$AT"\t"$CG"\t"$CT"\t"$GT"\t"$AA"\t"$CC"\t"$GG"\t"$TT  >> Daniel/Results/polyN_patterns.txt
done