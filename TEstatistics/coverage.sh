#!/bin/bash
if [ $# -ne 6 ]
then
	echo -e "\nusage: $0 <rpmk.out> <features.bed> <min_ovlp> <genome_size(bp)> <feature_type cod> <chrom_prefix> \n"
	echo -e "Description: This script calculates the total length (coverage) and overlap for each TE family against different possible features types in BED format"
	echo -e "Input:
		rpmk.out:	RepeatMasker output.
		features.bed:	Genomic query feature coordinates in bed format.
		feature_type cod (integer):  	1 -> genes
			     	 		2 -> CDS
			     	 		3 -> Introns
		             	 		4 -> UTRs
		             	 		5 -> flanking
		             	 		6 -> intergenic"
	echo -e "software requirements: bedtools"
	exit
fi


case $5 in
	1)
	 f_type="genes"
	 ;;
	2)
	 f_type="CDS"
	 ;;
	3)
	 f_type="introns"
	 ;;
	4)
	 f_type="UTRs"
	 ;;
	5)
	 f_type="2kflank"
	 ;;
	6)
	 f_type="intergenic"
	 ;;
	*)
	 exit
	 ;;
esac

rpmk_out=$1
features_bed=$2
min_ovlp=$3
genome_size=$4
prefixChr=$6
function merge_ovlp_TEs_by_fam() { #merge any two overlaping TE instances if they belong to the same family
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


function merge_ovlp_features() { #merge any two overlaping features from the input
	sort -k1,1 -k2,2n $features_bed > $features_bed.sort
	bedtools merge -i $features_bed.sort -c 4,5,6 -o distinct,distinct,distinct > $features_bed.merge
	awk -v OFS="\t" '{print $1,$2,$3,$4,$5,"."}' $features_bed.merge |sort -k1,1 -k2,2n > $features_bed.merge.bed
	rm  $features_bed.merge $features_bed.sort
}

function intersect_features-repeats(){ 
	bedtools intersect -a $1 -b $2 -wao > all_intersect
	awk -v pfx="$prefixChr" '$7~pfx {print $0}' all_intersect | sort -nk 13,13 | awk -v var="$min_ovlp" '$13 >= var {print $0}' > TE_intersect  #only TEs whose intersection is >= 50 bp are kept
	rm all_intersect $1
}

function calculate_feature_coverage_by_family () {  #output format: TE_name | total_overlap_with_the_queried_features(bp)
	awk '{print $10,$13}' $1 > TEnames_plus_ovlp.tmp
	awk '{print $10}' $1 | sort | uniq > TEnames.tmp
	awk 'NR==FNR{Arr[$1]=Arr[$1]+$2} NR!=FNR {print $1,Arr[$1]}' TEnames_plus_ovlp.tmp TEnames.tmp >> TE_$f_type\_intersect
	rm $1 TEnames_plus_ovlp.tmp
}

function calculate_TE_family_size() {
grep -v "(.*)" all.merge.sort | grep -v "rich" | awk -v var="$min_ovlp" '($3-$2) >= var {print $4,$3-$2}' | sort -k 2,2n > $rpmk_out.bed.$min_ovlp
awk -v g_size="$genome_size" 'NR==FNR{Arr[$1]=Arr[$1]+$2} NR!=FNR {print $1,Arr[$1],(Arr[$1]/g_size)*100}' $rpmk_out.bed.$min_ovlp TEnames.tmp >> TE_family_size   #format: TE_name |  TE_coverage(bp) | TE_coverage (%)
rm all.merge.sort $rpmk_out.bed.$min_ovlp TEnames.tmp
}


#------------------------------ MAIN -------------------------------------------_#
##################################################################################



merge_ovlp_TEs_by_fam $rpmk_out   					#  ->  all.merge.sort

merge_ovlp_features  							#  ->  $features_bed.merge.bed 

intersect_features-repeats $features_bed.merge.bed all.merge.sort  	# -> TE_intersect

calculate_feature_coverage_by_family TE_intersect  			# -> TE_$f_type\_intersect

calculate_TE_family_size   						# -> TE_family_size
							#Calculate TE-family size (i.e. genomic coverage)

#make a summary table
awk -v OFS="\t" -v feature=$f_type 'BEGIN{print "TE_name","TE_coverage_bp","TE_coverage_%","Overlap_with_"feature"_bp","Overlap_with_"feature"_%","Complement_Overlap_bp"} NR==FNR{Arr[$1]=$2} NR!=FNR {print $1,$2,$3,Arr[$1],(Arr[$1]*100)/$2,$2-Arr[$1]}' TE_$f_type\_intersect TE_family_size > TE_$f_type\_intersect.tab

rm TE_$f_type\_intersect TE_family_size
