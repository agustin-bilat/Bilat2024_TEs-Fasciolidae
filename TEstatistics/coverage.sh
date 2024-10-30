#!/bin/bash

# Check for the correct number of arguments
if [ $# -ne 4 ]; then
    echo -e "\nUsage: $0 <rmask.out> <genes.bed> <genome_assembly_size(bp)> <chromosome_prefix>\n"
    echo -e "Description: This script calculates the total length (coverage) for each TE family and the fraction that overlaps with both genic and intergenic regions.\n"
    echo -e "Inputs:"
    echo -e "  rmask.out:            RepeatMasker output."
    echo -e "  genes.bed:           Gene coordinates in BED format."
    echo -e "  genome_assembly_size: Total length of the assembly in base pairs."
    echo -e "  chromosome_prefix:    The name prefix used for scaffolds, contigs, or chromosomes. Both 'rmask.out' and 'genes.bed' files must use the same names for these sequences."
    echo -e "\nSoftware Requirements: bedtools"
    exit
fi

# Assign input arguments to variables
rpmk_out=$1
features_bed=$2
genome_size=$3
prefixChr=$4

# Convert RepeatMasker to BED and merge any two overlaping TE instances if they belong to the same family
function merge_ovlp_TEs_by_fam() { 
	mkdir dir_temporal
	cd dir_temporal
	grep -v "Simple_repeat" ./../$1 | awk '$10 !~/rich/{print $0}' | awk '$11!~/RNA/{print}' |  awk '$10 !~/polypurine/ && NR>3 {print $5,$6-1,$7,$10,".",$9}' |awk -v OFS="\t" '$6 ~ /C/{print $1,$2,$3,$4,".","-"} $6 ~/\+/{print $1,$2,$3,$4,".",$6}' | sort -k1,1 -k2,2n > $1.bed  #convert rpmk.out to bed format and remove repeats others than TEs
	awk '{print $4}' $1.bed | sort |uniq > names.tmp
	while read p; do awk -v var="$p" '$4==var {print $0}' $1.bed > $p.fam.bed ; done < names.tmp #split file by TE name
	for i in [a-zA-Z]*.fam.bed ; do bedtools merge -i $i -c 4,5,6 -o distinct,distinct,distinct > $i.merge; done # bedtools merge
	cat *bed.merge > all.merge  #re-join all TEs families in a bedfile format 
	sort -k1,1 -k2,2n all.merge > all.merge.sort
	rm names.tmp *.merge [a-zA-Z]*.fam.bed $1.bed
	mv all.merge.sort ./../
	cd ..
	rmdir dir_temporal
}

# Merge overlaping features (genes)
function merge_ovlp_features() { 
	sort -k1,1 -k2,2n $features_bed > $features_bed.sort
	bedtools merge -i $features_bed.sort -c 4,5,6 -o distinct,distinct,distinct > $features_bed.merge
	awk -v OFS="\t" '{print $1,$2,$3,$4,$5,"."}' $features_bed.merge |sort -k1,1 -k2,2n > $features_bed.merge.bed
	rm  $features_bed.merge $features_bed.sort
}

# Intersect features and repeats
function intersect_features-repeats(){ 
	bedtools intersect -a $1 -b $2 -wao > all_intersect
	awk -v pfx="$prefixChr" '$7~pfx {print $0}' all_intersect | sort -nk 13,13 | awk '$13 >= 1 {print $0}' > TE_intersect 
	rm all_intersect $1
}

# Calculate TE family coverage overlapping genes
function calculate_feature_coverage_by_family () {
	awk '{print $10,$13}' $1 > TEnames_plus_ovlp.tmp
	awk '{print $10}' $1 | sort | uniq > TEnames.tmp
	awk 'NR==FNR{Arr[$1]=Arr[$1]+$2} NR!=FNR {print $1,Arr[$1]}' TEnames_plus_ovlp.tmp TEnames.tmp >> TE_genes_intersect
	rm $1 TEnames_plus_ovlp.tmp
}

# Calculate total TE family size
function calculate_TE_family_size() {
grep -v "(.*)" all.merge.sort | grep -v "rich" | awk '($3-$2) >= 1 {print $4,$3-$2}' | sort -k 2,2n > $rpmk_out.bed.tmp
awk -v g_size="$genome_size" 'NR==FNR{Arr[$1]=Arr[$1]+$2} NR!=FNR {print $1,Arr[$1],(Arr[$1]/g_size)*100}' $rpmk_out.bed.tmp TEnames.tmp >> TE_family_size   #format: TE_name |  TE_coverage(bp) | TE_coverage (%)
rm all.merge.sort $rpmk_out.bed.tmp TEnames.tmp
}


#------------------------------ MAIN -------------------------------------------_#
##################################################################################



merge_ovlp_TEs_by_fam $rpmk_out   					#  ->  all.merge.sort

merge_ovlp_features  							#  ->  $features_bed.merge.bed 

intersect_features-repeats $features_bed.merge.bed all.merge.sort  	# -> TE_intersect

calculate_feature_coverage_by_family TE_intersect  			# -> TE_genes_intersect

calculate_TE_family_size   						# -> TE_family_size
							#Calculate TE-family size (i.e. genomic coverage)

#make a summary table
awk -v OFS="\t" -v feature="genes" 'BEGIN{print "TE_name","TE_coverage_bp","TE_coverage_%","Overlap_with_"feature"_bp","Overlap_with_"feature"_%","Complement_Overlap_bp"} NR==FNR{Arr[$1]=$2} NR!=FNR {print $1,$2,$3,Arr[$1],(Arr[$1]*100)/$2,$2-Arr[$1]}' TE_genes_intersect TE_family_size > TE_genes_intersect.tab

rm TE_genes_intersect TE_family_size
