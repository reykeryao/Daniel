#!bin/bash

# Pipeline for PE reads

# Need JA number, and part of the sample name from the RAW fastq file as input parameters
# Such as JA15599, and YQ1
# $1 is the JA number start with JA, and $2 is the sample name, better to choose a short but not an ambiguous one
# $3 is for UMI processing, if you have UMI, put UMI as $3, and default UMI is 10nt, change it in the script if your UMI is different
# 1. check folder structure, if not exit, make one
# folder structure: 
# RAW fastq is stored under $HOME/NGS/Data/$1, which should exist before running this script
# The folder will be checked and created while not existing is $HOME/NGS/Work/$1
# $HOME/NGS/Work/$1-----Trim
#                 |-----$2 (your sample name)-----Pass1
#                 |                         |-----Pass2
#                 |                         |-----Pass3
#                 |                         |-----Combined
#                 |-----Counts
# no "--norc"

JA=$1
Sample=$2
UMI=${3:-noUMI}

RAWfolder="$SCRATCH/Data/$JA"
Workfolder="$SCRATCH/Work/$JA"
Trimfolder="$SCRATCH/Work/$JA/Trim"
Samplefolder="$SCRATCH/Work/$JA/$Sample"
Countsfolder="$SCRATCH/Work/$JA/Counts"
if [ ! -d "$RAWfolder" ] ; then echo "This folder doesn't exist."; exit 1; fi
if [ ! -d "$Workfolder" ] ; then
        mkdir -p "$Trimfolder"
        mkdir -p "$Samplefolder/Pass1"
        mkdir "$Samplefolder/Pass2"
        mkdir "$Samplefolder/Combined"
        mkdir "$Countsfolder"
else
        if [ ! -d "$Samplefolder" ] ; then
                mkdir -p "$Samplefolder/Pass1"
                mkdir "$Samplefolder/Pass2"
                mkdir "$Samplefolder/Combined"
 else
                if [ ! -d "$Samplefolder/Pass1" ] ; then
                        mkdir "$Samplefolder/Pass1"
                fi
                if [ ! -d "$Samplefolder/Pass2" ] ; then
                        mkdir "$Samplefolder/Pass2"
                fi
                if [ ! -d "$Samplefolder/Combined" ] ; then
                        mkdir "$Samplefolder/Combined"
                fi
        fi
        if [ ! -d "$Trimfolder" ] ; then
                mkdir "$Trimfolder"
        fi
        if [ ! -d "$Countsfolder" ] ; then
                mkdir "$Countsfolder"
        fi
fi
#2 Start adapt trimming, UMI processing, default is no UMI
file="$RAWfolder/$Sample*.gz"
file1=$(echo $file|cut -f 1 -d " ")
file2=$(echo $file|cut -f 2 -d " ")
num_of_file=$(ls $file|wc -l)
if [ $num_of_file != 2 ];then echo "More than 2 files share the same sample name, make sure the name is uniq and you have PE reads"; exit 1;fi
trimed1="$Trimfolder/$Sample.1.fq.gz"
trimed2="$Trimfolder/$Sample.2.fq.gz"
if [ $UMI != "UMI" ] ; then
#        trimed_file2="$RAWfolder/$Sample.2.temp.fq.gz"
#        gunzip -c $file2 | fastx_trimmer -z -f 2 -i - -o $trimed_file2
        cutadapt -m 15 -O 3 --nextseq-trim=20 -j 48 --trim-n -q 20 -a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC  -g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGA -o $trimed1 -p $trimed2 $file1 $file2 1>"$Samplefolder/cutadapt.log"
        rm $RAWfolder/$Sample.2.temp.fq.gz
else
        dedupe_file1="$RAWfolder/$Sample.1.dedupe.fq.gz"
        dedupe_file2="$RAWfolder/$Sample.2.dedupe.fq.gz"
        clumpify.sh dedupe subs=1 in1=$file1 in2=$file2 out1=$dedupe_file1 out2=$dedupe_file2 2>"$Samplefolder/dedupe.log"
        trimed_file1="$RAWfolder/$Sample.1.temp.fq.gz"
        trimed_file2="$RAWfolder/$Sample.2.temp.fq.gz"
        gunzip -c $dedupe_file1 | fastx_trimmer -z -f 11 -i - -o $trimed_file1
        gunzip -c $dedupe_file2 | fastx_trimmer -z -f 2 -i - -o $trimed_file2
        trimed_file3="$RAWfolder/$Sample.3.temp.fq.gz"
        cutadapt -m 15 -O 3 --nextseq-trim=20 -j 48 --trim-n -q 20 -a AAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -g GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTT -A GATCGTCGGACTGTAGAACTCTGAACGTGTAGA -o $trimed1 -p $trimed_file3 $trimed_file1 $trimed_file2 1>"$Samplefolder/cutadapt.log"
        gunzip -c $trimed_file3 | fastx_trimmer -z -t 8 -i - -o $trimed2
        rm $RAWfolder/$Sample.*.temp.fq.gz
fi

#3 Pass1: hisat2 mapping
cd "$Samplefolder/Pass1"
rm -f *.*
hisat2 -p 48 -k 1 --rdg 1,3 --rfg 1,3 --mp 2,1 --no-mixed --no-discordant --no-spliced-alignment --new-summary --seed 740714 -x $LCREF/Ecoli_SC366/SC366_hisat -1 $trimed1 -2 $trimed2 2>Pass1.log |samtools view -bS -> Pass1.bam

samtools view -f4 Pass1.bam > unmapped.sam
samtools view -H Pass1.bam > filter.sam
samtools view -F4 Pass1.bam | bash $WORK/pipeline/072820_simple_filter
samtools view -bS filter.sam > filter.bam 
samtools fastq -N -1 unmapped.1.fq.gz -2 unmapped.2.fq.gz unmapped.sam
#samtools view -bS unmapped.sam | bam2fastx -APQN -o unmapped.fq.gz -
rm *.sam

#4 Pass2: bowtie2 mapping
cd "$Samplefolder/Pass2"
rm -f *.*
bowtie2 -p 48 -k 1 --rdg 1,3 --rfg 1,3 --mp 4 --ma 1 --no-mixed --no-discordant --very-sensitive-local --seed 740714 -x $LCREF/Ecoli_SC366/SC366_bowtie -1 ../Pass1/unmapped.1.fq.gz -2 ../Pass1/unmapped.2.fq.gz 2>Pass2.log | samtools view -bS - >Pass2.bam

samtools view -f4 Pass2.bam > unmapped.sam
samtools view -H Pass2.bam > filter.sam
samtools view -F4 Pass2.bam | bash $WORK/pipeline/072820_simple_filter
samtools view -bS filter.sam > filter.bam
samtools fastq -N -1 unmapped.1.fq.gz -2 unmapped.2.fq.gz unmapped.sam
#samtools view -bS unmapped.sam | bam2fastx -APQN -o unmapped.fq.gz -
rm *.sam

#5 Combined reads
cd "$Samplefolder/Combined/"
rm -f *.*
samtools merge primary.bam ../Pass1/filter.bam ../Pass2/filter.bam
samtools sort primary.bam > $Sample.sort.bam
samtools index $Sample.sort.bam

samtools view -h $Sample.sort.bam | awk '{if (($1~/^@/) || ($2==99) || ($2==147) || ($2==355) || ($2==403)) print}' | samtools view -bS - > $Sample.plus.sort.bam
samtools index $Sample.plus.sort.bam

samtools view -h $Sample.sort.bam | awk '{if (($1~/^@/) || ($2==83) || ($2==163) || ($2==339) || ($2==419)) print}' | samtools view -bS - > $Sample.minus.sort.bam
samtools index $Sample.minus.sort.bam

#6 Get the counts
bedtools bamtobed -mate1 -bedpe -i primary.bam | awk '{FS="\t"; OFS="\t"; if ($2>$5) $2=$5; if ($3<$6) $3=$6; print $1,$2,$3,$7,0,$9}' > primary.bed
bedtools coverage -s -counts -F 0.1 -a $LCREF/Ecoli_SC366/SC366.bed -b primary.bed  > "$Countsfolder/$Sample.counts"
