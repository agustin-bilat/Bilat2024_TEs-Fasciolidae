The manual curation workflow is based on the guidelines indicated at: Goubert, C., Craig, R.J., Bilat, A.F. et al. A beginner’s guide to manual curation of transposable elements. Mobile DNA 13, 7 (2022). https://doi.org/10.1186/s13100-021-00259-7.

The main steps and commands are indicated below:

### *De novo* search of TEs from genomes' assemblies ###
**Input:**  
Genomes' assemblies (FASTA) of *Fa. hepatica*, *Fa. gigantica* and *Fp. buski*.  
**Output:**  
RM2-libraries (set of repeats' consensus sequences in MULTI-FASTA format)  
**Software:**  
Raassssssssssssssssssssssssssssssssssssssssss


- RepeatModeler (RM2) (v 2.0.2)  
  `BuildDatabase -name <genome_assembly_name>.db <genome_assembly_name>.fna`
    
  `RepeatModeler -pa <number_of_parallel_search_jobs> -LTRStruct -database <genome_assembly_name>.db 2>rm2.err 1>rm2.out`

  
### Manual Curation of the *de novo* identified repeats ###

**Input:**  
A RM2-library (MULTI-FASTA)

**Output:**  
A curated TE libraries (MULTI-FASTA)

1. `cd-hit-est -i <RM2-library_name>.fa -o <Curated library_name>.fa.cdhit -c 0.8 -n 5 -G 0 -aS 0.8 -d 0 -g 1 -b 500`

2. `faSplit byname <library_name>.fa.cdhit <Directory_name>/`

Within each directory named as <Directory_name> we get the genomic insertions of each family by using the script [mkfasta_from_megablast.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/mkfasta_from_megablast.sh):

3. a) `for i in *fam*.fa; do bash mkfasta_from_megablast.sh <genome_assembly_name>.fna $i 0 1000 1000 <chromsizes>.fa <genome_assembly_database_prefix> ; done`

For the families in which less than 50 blast hits were obtained, the script [mkfasta_fromBlast.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/mkfasta_from_blastn.sh) is run as follows:

3. b) `for i in *fam*.fa; do bash mkfasta_fromBlast.sh <genome_assembly_name>.fna $i 0 1000 1000 <chromsizes>.fa <genome_assembly_database_prefix> ; done`

The set of families having at least 50 hits (*bed.fa files) are moved into a new directory were the script [ready_for_MSA.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/ready_for_MSA.sh) is run to subset a maximum number of sequences for each family to be aligned.

4. `for i in *fam*bed.fa ; do bash /home/agustin/Programs/TE_scripts_ab2112/TEcuration/ready_for_MSA.sh $i 200 50 ; done`

Multiple sequence alignment with [MAFFT](https://mafft.cbrc.jp/alignment/software/):

5. `for i in *.fa ; do mafft --thread 4 $i > $i.maf ; done`

The alignments are manually edited with _aliview_ as described by Goubert, C., 2022. For families were the borders of the alignment are not reached, we repeat the steps 3 to 5 changing the values of the command in step 3 in order to increase the flanking sequence of the blast hits to be extracted as fasta (in the example 1 kb at either side).  

When manual edition is ready, further automatic edition is made with the _CIAlign_ software with the following options and parameters:  

6. `for i in *.maf; do CIAlign --infile $i --remove_insertions --insertion_min_flank 3 --insertion_min_size 1 --insertion_max_size 500 --remove_divergent --remove_divergent_minperc 0.65 --crop_ends --remove_short --remove_min_length 100 --plot_input --plot_output --outfile_stem $i ; done`

We check the mini plots of the curated alignments for further manual curation when necessary (for example, if insertions are still observed in an aligmnet we eliminate them with aliview). Additional sequences might be eliminated based on the observation of the aligmnets plots (higly-fragmented/unconserved and low-copy aligments after CIAlign trimming).  
  
Create consensus sequence from each curated alignment (EMBOSS:6.6.0.0):

7. `for i in <repeat_alignments> ; do cons -sequence $i -outseq $i.cons.fa -plurality 0.1 -name $i ; done`  

For each species, the consensus are merged in multi-fasta files and _cdhit_ is run as in the step one of manual curation to eliminate redundancies.  The representative consensus for each cluster are kept in the library. 

8. Structural Characterization of Consensus:
- [TE-Aid](https://doi.org/10.1186/s13100-021-00259-7)  
  `for i in *cons.fa; do bash TE-Aid --query $i --genome <genome_assembly_name>.fna -m 400 -e 10e-10 -o TEAid_out`
- Protein domains (EMBOSS:6.6.0.0 and pfam_scan.pl tools are required)  
  `transeq -sequence <species_cdhit>.fa -outseq <species_cdhit>.fa.transeq -frame 6 -clean`  
  `pfam_scan.pl -fasta <species_cdhit>.fa.transeq -dir <pfam_database_path> -outfile <species>_pfams.txt`  

A Neighbour-Joining (NJ) tree of aligned Reverse Transcriptase (RVT) domains from all the curated consensus was generated with MEGA11 (Tamura et al., 2021) to get an global overview of the types of retrotranposons.  
Based on both the NJ-tree and the structural characterization, TEs were classified as DNA transposons (class II) or as LINE, LTR or PLE (class I).

9. Phylogenetic characterization with iqtree (v 1.6.12):

Align the RVT or DDE/Tase identified within translated consensus:
`mafft --localpair --maxiterate 1000 <protein-domain.aa> > <protein-domain.aa>.maf`
Make maximum-likelihood trees:
`iqtree -s <protein-domain_seq>.maf -m MFP -nt 2 -bb 1000 -bnni`  

Sequences are renamed and classified to finally obtain the curated libraries (see article).


aa
