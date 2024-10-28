#!/bin/bash
if [ $# -ne 1 ]
then
        echo -e "\nusage: $0 <out_dir_name>\n"
        echo "move and rename files from the output directory after runing insertions_to_fasta.sh into a single directory.\n This should be run within the directory in which the former script was run"
        exit
fi

out_dir=$1


for d in */; do
cd $d
cp total_insertions.fasta ${d%/}.insertions.fas
mv ${d%/}.insertions.fas ..
cd ..
done

mkdir $out_dir
mv *insertions.fas $out_dir/
