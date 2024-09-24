#!/bin/bash

if [ $# -ne 7 ]
then
    echo -e "\nusage: $0 <genome> <fasta.in> <min_length> <flank_l> <flank_r> <chrom_sizes> <genome_fasta> \n"
    echo -e "DESCRIPTION: This script takes a fasta sequence <fasta.in>, blasts it to the genome, recovers locations with alingment length > <min_length>, prints them as bed file, extends bed coordinates in each direction, and makes fasta from that BED.\n"

    echo -e "INPUT:       <fasta.in>    query sequence if fasta format"
    echo -e "             <min_length>  min length of the blast hit. If set to 0, min length = (length of query)/2"
    echo -e "             <flank>       number of bases to extend the genome coordinates of the matched locus upstream (flank_l) and downstream (flank_r) relative to the plus "+" strand. If only one direction is needed to be extended, just set to 0 the other flank (the one you don't want to extend)." 
    
    echo -e "OUTPUT:      produces a <fasta.in.bed> which is the blast results BED file file"
    echo -e "             produces a <fasta.in.blast.flank.bed> which is the extended BED"
    echo -e "             produces a <fasta.in.blast.flank.bed.fa> bedtools getfasta \n"
     
    exit
fi

genome=$1
fasta_in=$2
out=`basename $2`
min_length=$3
flank_l=$4
flank_r=$5
chrom_sizes=$6
genome_fasta=$7


if [ $3 == 0 ]
then 
	min_length=`grep -v ">" $2 | wc | awk '{print int($3/2)}'`
else
	min_length=$3
fi

echo "query sequence" $out
echo "minimum length of blast hit = " $min_length
echo "the hit locus will be extended" = $flank "bases in each direction"

echo "#qseqid sseqid pident length mismatch qstart qend sstart send sstrand" > $out.blast.o
blastn -query $fasta_in -db $genome -outfmt "6 qseqid sseqid pident length mismatch qstart qend sstart send sstrand" -evalue 1e-20 | awk -v "ml=$min_length" '{OFS="\t"; if ($4 > ml) {print $0}}' >> $out.blast.o

awk '{OFS="\t"; if ($1~/^\#/) {} else { if ($10~/plus/) {print $2, $8, $9, $1, $3, "+"} else {print $2, $9, $8, $1, $3, "-"}}}' < $out.blast.o > $out.blast.bed

bedtools slop -s  -i $out.blast.bed  -g $chrom_sizes -l $flank_l -r $flank_r > $out.blast.flank.bed

bedtools getfasta -fi  $genome_fasta -fo $out.blast.flank.bed.fa  -bed $out.blast.flank.bed -s 

fasta_count=`grep -c ">" $out.blast.flank.bed.fa`

echo "the fasta has "$fasta_count " sequences"

#rm $out.blast.o $out.blast.bed 
