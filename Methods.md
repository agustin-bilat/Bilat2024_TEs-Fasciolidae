1. ### *De novo* search of TEs from genomes' assemblies ###

- RepeatModeler (v 2.0.2)\
  `BuildDatabase -name <genome_assembly_name>.db <genome_assembly_name>.fna`
    
  `RepeatModeler -pa <number_of_parallel_search_jobs> -LTRStruct -database <genome_assembly_name>.db 2>rm2.err 1>rm2.out`
  
- SINE_scan

2. ### Manual Curation of the *de novo* identified repeats ###
