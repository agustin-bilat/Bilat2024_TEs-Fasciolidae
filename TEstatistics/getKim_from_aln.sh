#!/bin/bash

# Check arguments
if [ $# -ne 4 ]; then
    echo -e "\nUsage: $0 <chromosome/scaffold_prefix_name> <alignment (.aln) file> <total assembly size (bp)> <assembly_name>"
    exit
fi

# Assign input arguments to variables
prefix=$1              # Prefix for chromosome/scaffold
align=$2               # Alignment file
genome_size=$3        # Total assembly size in base pairs
assembly=$4            # Assembly name

# Filter alignment file and format it for further processing
grep -v "Simple_repeat" $align |  grep -E "$prefix.*#|Kimura.*CpG" | awk -v OFS="\t" -v pfix="$prefix" '($0 ~ pfix) && ($9!="C") \
{chrom=$5;start=($6-1);end=$7;rep_name=$9;sense="+"} ($0~ pfix) && ($9=="C"){chrom=$5;start=($6-1);end=$7;rep_name=$10;sense="-"} \ 
/Kimura/{print chrom,start,end,rep_name,".",sense,$5}' | grep -v "Satellite" > kimura.TE

# Create a file with TE name, size, and kimura-bins
awk -v OFS="\t" '{print $4, $3 - $2, $7}' kimura.TE | \
    awk -F '.' '{print $1}' | \
    awk -v OFS="\t" '{print $1, $2, $3}' > rpmask_kimura-and-size  

# Split output by kimura bin
awk -v out_file="bin_kimura.out" '{print $0 >> $3"."out_file; close($3"."out_file)}' rpmask_kimura-and-size

# Generate unique names for each output file
for i in *.out; do 
    awk '{print $1}' "$i" | sort | uniq > "$i.names"
done

# Process each kimura output file and generate results
for i in *kimura.out; do  
    awk 'NR==FNR {l[$1] = l[$1] + $2; k = $3} NR != FNR {print $1, l[$1], k}' "$i" "$i.names" > "$i.result"
done

# Combine results and sort them
cat *.out.result | sort -k1,1 -k3,3b > kimura.total_results

# Format the final results and output to a tab-separated file
awk -F '#' '{print $1, $2}' kimura.total_results | \
    sort -k1,1 -k4,4n | \
    awk -v OFS="\t" -v var1="$assembly" -v var2="$genome_size" '
        BEGIN { print "species", "TEname", "TEtype", "percent_of_genome", "kimura" } 
        { print var1, $1, $2, $3 * 100 / var2, $4 }' | \
    sed 's/,/./g' > "${assembly}_kim.tab"

# Clean up temporary files
rm *kimura.out* rpmask_kimura-and-size kimura.total_results kimura.TE
