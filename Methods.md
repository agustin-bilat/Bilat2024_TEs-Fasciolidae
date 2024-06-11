### *De novo* search of TEs from genomes' assemblies ###
**Input:**  
Genomes' assemblies (FASTA) of *Fa. hepatica*, *Fa. gigantica* and *Fp. buski*.  
**Output:**  
rm2 libraries libraries (set of repeat consensus sequences in FASTA format)

- RepeatModeler (v 2.0.2)  
  `BuildDatabase -name <genome_assembly_name>.db <genome_assembly_name>.fna`
    
  `RepeatModeler -pa <number_of_parallel_search_jobs> -LTRStruct -database <genome_assembly_name>.db 2>rm2.err 1>rm2.out`

- SINE_scan (v 1.1.1)

  
### Manual Curation of the *de novo* identified repeats ###

**Input:**  
raw TE libraries of *Fa. hepatica*, *Fa. gigantica* and *Fp.buski*:

**Output:**  
(curated) TE libraries (set of curated TE consensus sequences in FASTA format)

1. `cd-hit-est -i <library_name>.fa -o <library_name>.fa.cdhit -c 0.8 -n 5 -G 0 -aS 0.8 -d 0 -g 1 -b 500`

2. `faSplit byname <library_name>.fa.cdhit <Directory_name>/`

Within each directory named as <Directory_name> we get the genomic insertions of each family by using the script [mkfasta_from_megablast.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/mkfasta_from_megablast.sh):

3. `for i in *fam*.fa; do bash mkfasta_fromBlast.sh <genome_assembly_name>.fna $i 0 1000 1000 <chromsizes>.fa <genome_assembly_database_prefix> ; done`

For the families in which less than fifty blast hits were obtained, the script [mkfasta_fromBlast.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/mkfasta_from_blastn.sh) is run as follows:

4. `for i in *fam*.fa; do bash mkfasta_fromBlast.sh <genome_assembly_name>.fna $i 0 1000 1000 <chromsizes>.fa <genome_assembly_database_prefix> ; done`

The set of families having at least 50 hits (*bed.fa files) are moved into a new directory were the following script is run [ready_for_MSA.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/ready_for_MSA.sh): 

5. `for i in *fam*bed.fa ; do bash /home/agustin/Programs/TE_scripts_ab2112/TEcuration/ready_for_MSA.sh $i 200 50 ; done`

Multiple sequence alignment with [MAFFT](https://mafft.cbrc.jp/alignment/software/):

6. `for i in *.fa ; do mafft --thread 4 $i > $i.maf ; done`

After making some manual edition of the alignments with _aliview_ (see the paper's methods section) further automatic edition is made with the _CIAlign_ software with the following options and parameters: 

7. `for i in *.maf; do CIAlign --infile $i --remove_insertions --insertion_min_flank 3 --insertion_min_size 1 --insertion_max_size 500 --remove_divergent --remove_divergent_minperc 0.65 --crop_ends --remove_short --remove_min_length 100 --plot_input --plot_output --outfile_stem $i ; done`

Create consensus sequence from each curated alignment (EMBOSS:6.6.0.0):

8. `for i in <repeat_alignments> ; do cons -sequence $i -outseq $i.cons.fa -plurality 0.1 -name $i ; done`  

For each species, the consensus are merged in multi-fasta files and _cdhit_ is run as in the step one of manual curation. The representative consensus for each cluster are kept for the following analysis. 

9. Structural Characterization of Consensus:
- [TE-Aid](https://doi.org/10.1186/s13100-021-00259-7)  
  `for i in *cons.fa; do bash TE-Aid --query $i --genome <genome_assembly_name>.fna -m 400 -e 10e-10 -o TEAid_out`
- Protein domains (EMBOSS:6.6.0.0 and pfam_scan.pl tools are required)  
  `transeq -sequence <species_cdhit>.fa -outseq <species_cdhit>.fa.transeq -frame 6 -clean`  
  `pfam_scan.pl -fasta <species_cdhit>.fa.transeq -dir <pfam_database_path> -outfile <species>_pfams.txt` 
    
10. Phylogenetic characterization:


