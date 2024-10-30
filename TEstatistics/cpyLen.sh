#!/bin/bash
if [ $# -ne 2 ]
then
        echo -e "\nusage: $0 <repeat_masker.out> <assembly_name (short)> \n"
        echo -e "Description: Obtain the length of TE insertions (copy-length) from RepeatMasker output table for each TE family.\nOutput format: assembly_name\tTEname\tTEtype\tlength" 
        exit
fi

repMask=$1
assName=$2


# Convert RepeatMasker.out table into:  chromosome | chr-beg | chr-end | fam-name | "." | str | TE-type | insertion-beg | insertion-end | ID 
awk '$11!="Unknown"{print}' $repMask | grep -v "Simple_repeat" | grep -v "Low_complexity" | awk '$11!="Satellite"{print}' | awk -v OFS="\t" 'NR>3 && $9~/C/{print $5,$6-1,$7,$10,".","-",$11,$14-1,$13,$15}NR>3 && $9~/+/{print $5,$6-1,$7,$10,".","+",$11,$12-1,$13,$15}' |  sort -k10,10n -k8,8n > table1.tmp

# Get length of TE copies ( TEname#ID | TEtype | copy-length)
awk -v OFS="\t" '{print $4"#"$10,$7,$3-$2}' table1.tmp > table2.tmp

# Get unique lines by (TEname AND ID):
awk '{print $1}' table2.tmp | sort | uniq > names.tmp

# Sum the length of TE copies by ID:
awk 'NR==FNR{len[$1]+=$3;type[$1]=$2}NR!=FNR{print $1,type[$1],len[$1]}' table2.tmp names.tmp > defrag_TEins.tab

# Reformat table and add assembly name:
awk -F '#' '{print $1,$2}' defrag_TEins.tab | awk -v OFS="\t" -v var="$assName" '{print var,$1,$3,$4}' > $2.defragment_TEins.tab
sort -k1,1 -k2,2 -k4,4nr $2.defragment_TEins.tab > $2.TEinsertions.tab


rm table1.tmp table2.tmp names.tmp defrag_TEins.tab $2.defragment_TEins.tab
