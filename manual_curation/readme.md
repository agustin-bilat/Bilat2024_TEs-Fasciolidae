This file describes the specific steps and parameters associated with the manual curation of repeats that is mentioned in Methods section at ***the article***. The metodology is based on the guidelines described elsewhere (Goubert, C., Craig, R.J., Bilat, A.F. et al. A beginner’s guide to manual curation of transposable elements. Mobile DNA 13, 7 (2022). https://doi.org/10.1186/s13100-021-00259-7). The bibliography associated to each software can be found in that article. The scripts indicated below are included in modified from the ones present in the Goubert article.

## Abbreviations ##  

***RM2-library:*** RepeatModeler (v 2.0.2) output of multi-fasta consensus sequences.  

  
## Manual Curation ##

**Input:**  
RM2-library (MULTI-FASTA)

**Output:**  
curated TE-library (MULTI-FASTA)

1. `cd-hit-est -i <RM2-library_name>.fa -o <RM2-library_name>.fa.cdhit -c 0.8 -n 5 -G 0 -aS 0.8 -d 0 -g 1 -b 500`

2. `faSplit byname <RM2-library_name>.fa.cdhit <Directory_name>/`

Within the directory "<Directory_name>" the script [mkfasta_from_megablast.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/mkfasta_from_megablast.sh) is run as follows:

3. a) `for i in *fam*.fa; do bash mkfasta_from_megablast.sh <genome_assembly_name>.fna $i 0 1000 1000 <chromsizes>.fa <genome_assembly_database_prefix> ; done`

Output families (*bed.fa, MULTI-FASTA files) generated with at least 50 megablast hits are moved into a new directory until step 4.  
For the remaining families the script output files are removed, and the original RM2-queries are reused as inputs into a new more flexible script ([mkfasta_fromBlast.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/mkfasta_from_blastn.sh)) as it is described below:

3. b) `for i in *fam*.fa; do bash mkfasta_fromBlast.sh <genome_assembly_name>.fna $i 0 1000 1000 <chromsizes>.fa <genome_assembly_database_prefix> ; done`

Output families (*bed.fa, MULTI-FASTA files) with at least blastn hits are moved into the same directory were the other high copy megablast hits families were located.  
Next, the script [ready_for_MSA.sh](https://github.com/agustin-bilat/Bilat2024_TEs-Fasciolidae/blob/main/scripts/ready_for_MSA.sh) is run within that directory to subset a maximum number of 200 sequences:

4. `for i in *fam*bed.fa ; do bash ready_for_MSA.sh $i 200 50 ; done`

Multiple sequence alignment are then generated using [MAFFT](https://mafft.cbrc.jp/alignment/software/):

5. `for i in *.fa ; do mafft --thread 4 $i > $i.maf ; done`

The alignments are manually edited with _aliview_ as described by Goubert, C., 2022. For families were the borders of the alignment are not reached, the steps 3 to 5 are repeated by changing the parameters in step 3 in order to increase the flanking sequences of hits to be extracted as fastas with respect to the first iteration (in which 1000 bp at either side of each hit was extended).  

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