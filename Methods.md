1. ### *De novo* search of TEs from genomes' assemblies ###
**Input:**  
Genomes' assemblies (FASTA) of *Fa. hepatica*, *Fa. gigantica* and *Fp. buski*.  
**Output:**  
rm2 libraries libraries (set of repeat consensus sequences in FASTA format)

- RepeatModeler (v 2.0.2)  
  `BuildDatabase -name <genome_assembly_name>.db <genome_assembly_name>.fna`
    
  `RepeatModeler -pa <number_of_parallel_search_jobs> -LTRStruct -database <genome_assembly_name>.db 2>rm2.err 1>rm2.out`

- SINE_scan (v 1.1.1)

  
2. ### Manual Curation of the *de novo* identified repeats ###

**Input:**  
raw TE libraries of *Fa. hepatica*, *Fa. gigantica* and *Fp.buski*:

**Output:**  
(curated) TE libraries (set of curated TE consensus sequences in FASTA format)

1. `cd-hit-est -i <library_name>.fa -o <library_name>.fa.cdhit -c 0.8 -n 5 -G 0 -aS 0.8 -d 0 -g 1 -b 500`

2. `faSplit byname <library_name>.fa.cdhit <Directory_name>/`

Within each directory named as <Directory_name> we get the genomic insertions of each family by using the script [mkfasta_from_megablast.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/mkfasta_from_megablast.sh):

3. `for i in *fam*.fa; do bash mkfasta_fromBlast.sh <genome_assembly_name>.fna $i 0 1000 1000 <chromsizes>.fa <genome_assembly_database_prefix> ; done`

For the families in which less than fifty blast hits were obtained, the following script is run [mkfasta_fromBlast.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/mkfasta_from_blastn.sh)

4. `for i in *fam*.fa; do bash mkfasta_fromBlast.sh <genome_assembly_name>.fna $i 0 1000 1000 <chromsizes>.fa <genome_assembly_database_prefix> ; done`

The set of consensus having at least 50 hits are located in a new directory (files with the suffix name "*bed.fa") and the following script is run [ready_for_MSA.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/ready_for_MSA.sh) 

5. `for i in *fam*bed.fa ; do bash /home/agustin/Programs/TE_scripts_ab2112/TEcuration/ready_for_MSA.sh $i 200 50 ; done`

Multiple sequence alignment with [MAFFT](https://mafft.cbrc.jp/alignment/software/)

6. `for i in *.fa ; do mafft --thread 4 $i > $i.maf ; done 
