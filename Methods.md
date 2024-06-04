1. ### *De novo* search of TEs from genomes' assemblies ###
**Input:**  
Genomes' assemblies (FASTA) of *Fa. hepatica*, *Fa. gigantica* and *Fp. buski*.  
**Output:**  
raw TE libraries (set of repeat consensus sequences in FASTA format)

- RepeatModeler (v 2.0.2)  
  `BuildDatabase -name <genome_assembly_name>.db <genome_assembly_name>.fna`
    
  `RepeatModeler -pa <number_of_parallel_search_jobs> -LTRStruct -database <genome_assembly_name>.db 2>rm2.err 1>rm2.out`

- SINE_scan

2. ### Manual Curation of the *de novo* identified repeats ###

**Input:**  
raw TE libraries *Fa. hepatica*, *Fa. gigantica* and *Fp.buski*:

**Output:**  
(curated) TE libraries (set of curated TE consensus sequences in FASTA format)

1. From the initial
