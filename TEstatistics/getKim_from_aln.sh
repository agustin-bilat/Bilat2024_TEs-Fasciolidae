#!/bin/bash

if [ $# -ne 4 ]
then
	echo -e "\nusage: $0 <chromosome/scaffold_prefix_name> <alignment (.aln) file> <total assembly size (bp)> <assembly_name>"

	exit
fi

prefix=$1
align=$2
genome_size=$3
assembly=$4

grep -v "Simple_repeat" $align |  grep -E "$prefix.*#|Kimura.*CpG" | awk -v OFS="\t" -v pfix="$prefix" '($0 ~ pfix) && ($9!="C") {chrom=$5;start=($6-1);end=$7;rep_name=$9;sense="+"} ($0~ pfix) && ($9=="C"){chrom=$5;start=($6-1);end=$7;rep_name=$10;sense="-"} /Kimura/{print chrom,start,end,rep_name,".",sense,$5}' | grep -v "Satellite" > kimura.TE

awk -v OFS="\t" '{print $4,$3-$2,$7}' kimura.TE | awk -F '.' '{print $1}' | awk -v OFS="\t" '{print $1,$2,$3}' > rpmask_kimura-and-size	# name | size | kimura-bins

awk -v out_file="bin_kimura.out" '{print $0 >> $3"."out_file; close($3"."out_file)}' rpmask_kimura-and-size

for i in *.out;do awk '{print $1}' $i | sort | uniq > $i.names ; done

for i in *kimura.out; do  awk 'NR==FNR{l[$1]=l[$1]+$2;k=$3}NR!=FNR{print $1,l[$1],k}' $i $i.names > $i.result ; done

cat *.out.result | sort -k1,1 -k3,3b > kimura.total_results

awk -F '#' '{print $1,$2}' kimura.total_results | sort -k1,1 -k4,4n | awk -v OFS="\t" -v var1="$assembly" -v var2="$genome_size" 'BEGIN{print "species","TEname","TEtype","percent_of_genome","kimura"}{print var1,$1,$2,$3*100/var2,$4}' | sed 's/,/./g' > ${assembly}_kim.tab

rm *kimura.out* rpmask_kimura-and-size kimura.total_results kimura.TE
