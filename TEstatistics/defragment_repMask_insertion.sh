#!/bin/bash
if [ $# -ne 3 ]
then
        echo -e "\nusage: $0 <rep.mask.out> <genome.fna> <families_length>\n"
        echo "Select and defragment 500 random insertion sequences for each TE family from the rep.masker .out file and then create the multiple sequence alignments for each family from those defragmented sequences. Only defragmented insertions larger than half the length of the family's consensus are taken for the alignment"
        exit
fi

repMask=$1
genome=$2
fams_len=$3


function get_fastas() {

#4a/ randomly select 500 insertions' IDs
awk '{print $10}' bedfile | sort | uniq | sort -R | head -500 > random_500.IDs
awk 'NR==FNR{id[$1]=$1}NR!=FNR && ($10 in id){print $0}' random_500.IDs bedfile > random_500.bed

#4b/ get IDs of insertions with multiple and single rep. masker hits within those 500:

awk '{print $10}' random_500.bed | sort | uniq -d > random_bed_multi.IDs
awk '{print $10}' random_500.bed | sort | uniq -u > random_bed_single.IDs

rm random_500.IDs random_500.bed

#4c/ get fastas for single-hit insertions:

# if file is not empty ; then do
single=random_bed_single.IDs
if [ -s "${single}" ] ; then
	awk -v OFS="\t" 'NR==FNR{arr[$1]=$1}NR!=FNR && ($10 in arr) {print $1,$2,$3,$4"#"$10,$5,$6 >> "single_insertions.Bed"}' random_bed_single.IDs bedfile

	sort -k1,1 -k2,2n single_insertions.Bed > single_insertions.Bed6

	bedtools getfasta -fi ./../$genome -fo total_insertions.fasta -bed single_insertions.Bed6 -s -nameOnly

	rm single_insertions.Bed single_insertions.Bed6 
fi


#4d/ get the defragmented fastas from multiple-hit insertions:
multi=random_bed_multi.IDs
if [ -s "${multi}" ]; then
	awk -v OFS="\t" 'NR==FNR{arr[$1]=$1}NR!=FNR && ($10 in arr) {print $1,$2,$3,$4"#"$10,$5,$6 >> $4"#"$10".Bed.tmp"}' random_bed_multi.IDs bedfile

	for i in *.Bed.tmp; do bedtools getfasta -fi ./../$genome -fo $i.fa -bed $i -s -nameOnly ;done
	for i in *.Bed.tmp.fa ;do union -sequence $i -outseq $i.union ; done

	cat *.fa.union >> total_insertions.fasta
	rm *.Bed.tmp*
fi

}





###################  MAIN #######################



#1/ from $repMask and $fam_len get a bed-like file of TEs formated as follows: 
# chromosome | chr-beg | chr-end | fam-name | "." | str | TE-type | insertion-beg | insertion-end | insertion-ID | Family(cons)-length

awk '$11!="Unknown"{print}' $repMask | grep -v "Simple_repeat" | grep -v "Low_complexity" | awk '$11!="Satellite"{print}' | awk -v OFS="\t" 'NR>3 && $9~/C/{print $5,$6-1,$7,$10,".","-",$11,$14-1,$13,$15}NR>3 && $9~/+/{print $5,$6-1,$7,$10,".","+",$11,$12-1,$13,$15}' |  sort -k10,10n -k8,8n > defragmented_TEs.bed

awk -v OFS="\t" 'NR==FNR{len[$1]=$2}NR!=FNR && ($4 in len){print $0,len[$4]}' $fams_len defragmented_TEs.bed > defragm_TEs.bed


#2/ remove insertions smaller than half the length of its family's consensus

awk 'NR==FNR{sum[$4"#"$10]+=$9-$8}NR!=FNR && sum[$4"#"$10]>=(int($11*0.5)){print}' defragm_TEs.bed defragm_TEs.bed > defragmTE_filtered.bed
rm defragm_TEs.bed defragmented_TEs.bed 


#3/ separate bed by family name in different directories

awk '{print $4}' defragmTE_filtered.bed | sort | uniq > names_fams.txt #fam names
while read p ;do mkdir $p ; done < names_fams.txt
while read p; do awk -v var="$p" '$4==var{print}' defragmTE_filtered.bed > $p/bedfile ; done < names_fams.txt

rm names_fams.txt


#4/ within each directory, get a multi-fasta file with at most 500 rdmly-selected defragmented sequences of the corresponding family.

for d in */; do
cd $d
get_fastas
cd ..
done
