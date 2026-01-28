for i in HX{01..04}
do
	mkdir -p Daniel/Merged/$i"_meme"
	cd Daniel/Merged/$i"_meme"
	zcat ../$i.fa.gz > $i.fa
	meme $i.fa -dna -oc . -nostatus -time 14400 -mod zoops -nmotifs 10 -minw 6 -maxw 50 -objfun classic -markov_order 0 -p 12
	cd ../../../
done
