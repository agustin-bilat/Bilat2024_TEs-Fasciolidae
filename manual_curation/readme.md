This readme file describes the specific steps and parameters used for the manual curation of repeats mentioned at the section Methods of the manuscript. The methodology is based on the guidelines described elsewhere (Goubert, C., Craig, R.J., Bilat, A.F. et al. A beginner’s guide to manual curation of transposable elements. Mobile DNA 13, 7 (2022). https://doi.org/10.1186/s13100-021-00259-7). The reference for each software can be found in manuscript. Some scripts that were used for manual curation were obtained by sligthly modifying the ones provided at the Goubert's et al. article and are included in this directory: [manual_curation/](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/manual_curation/).

## Abbreviations ##  

***RM2-library:*** RepeatModeler (v 2.0.2) output of multi-fasta consensus sequences.  

  
## Manual Curation ##

**Input:**  
RM2-library (MULTI-FASTA)

**Output:**  
curated TE-library (MULTI-FASTA)

1. `cd-hit-est -i <RM2-library_name>.fa -o <RM2-library_name>.fa.cdhit -c 0.8 -n 5 -G 0 -aS 0.8 -d 0 -g 1 -b 500`

2. `faSplit byname <RM2-library_name>.fa.cdhit <Directory_name>/`

Within the directory "<Directory_name>" the script [mkfasta_from_megablast.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/manual_curation/mkfasta_from_megablast.sh) is run as follows:

3. a) `for i in *fam*.fa; do bash mkfasta_from_megablast.sh <genome_assembly_name>.fna $i 0 1000 1000 <chromsizes>.fa <genome_assembly_database_prefix> ; done`

Outputs for families with at least 50 megablast hits are moved into a separate directory. For the remaining families, the script's output files are removed, and the original RM2-queries are re-used as inputs into a more flexible script (blastn instead of megablast) ([mkfasta_fromBlast.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/manual_curation/mkfasta_from_blastn.sh)) as it is described below:

3. b) `for i in *fam*.fa; do bash mkfasta_fromBlast.sh <genome_assembly_name>.fna $i 0 1000 1000 <chromsizes>.fa <genome_assembly_database_prefix> ; done`

Output files from families with at least 50 blastn hits are moved into the directory mentioned above, were the megablast output files should be located.  
Next, the script [ready_for_MSA.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/manual_curation/ready_for_MSA.sh) is run within that directory to subset a maximum number of 200 sequences:

4. `for i in *fam*bed.fa ; do bash ready_for_MSA.sh $i 200 50 ; done`

Multiple sequence alignment are then generated using [MAFFT](https://mafft.cbrc.jp/alignment/software/). For families having more than the maximum threshold number of selected sequences ("200" in the case shown in step 4) the alignment is made for the ready_for_MSA.sh's output *.fa files but not for the original multi-fasta files.

5. `for i in *.fa ; do mafft --thread 4 $i > $i.maf ; done`

The alignments are manually edited with _aliview_ as described by Goubert, C., 2022. For families were the borders of the alignment are not reached, the steps 3 to 5 are repeated by changing the parameters in step 3 in order to increase the flanking sequences of hits to be extracted as fastas with respect to the first iteration (in which 1000 bp at either side of each hit is extended). The corresponding older outputs are removed.

Next, further automatic edition is made with the _CIAlign_ software with the following options and parameters:  

6. `for i in *.maf; do CIAlign --infile $i --remove_insertions --insertion_min_flank 3 --insertion_min_size 1 --insertion_max_size 500 --remove_divergent --remove_divergent_minperc 0.65 --crop_ends --remove_short --remove_min_length 100 --plot_input --plot_output --outfile_stem $i ; done`

The input and output alignment's mini-plots are checked to evaluate if further manual curation is necessary (for example, if insertions are still observed in the aligmnet they are removed with aliview).   

Create consensus sequence from each curated alignment (EMBOSS:6.6.0.0):

7. `for i in <repeat_alignments> ; do cons -sequence $i -outseq $i.cons.fa -plurality 0.1 -name $i ; done`  

The collection of consensus are then merged in multi-fasta files and _cdhit_ is run again to eliminate redundancies. The representative consensus for each cluster are kept in the library for characterization, naming, classification and further filtering as it is described in the Methods section of the manuscript.
