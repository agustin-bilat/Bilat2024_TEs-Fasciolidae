1. ### *De novo* search of TEs from genomes' assemblies ###
**Input:**  
Genomes' assemblies (FASTA) of *Fa. hepatica*, *Fa. gigantica* and *Fp. buski*.  
**Output:**  
raw TE libraries (set of repeat consensus sequences in FASTA format)

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

Within the <Directory_name> directory, which should only contain the fasta sequences of one of the three libraries we run the following command:

3. `for i in \*fam\*.fa; do bash <4_mkfasta_fromBlast.sh>
