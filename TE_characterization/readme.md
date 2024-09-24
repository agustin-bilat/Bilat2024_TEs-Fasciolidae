1. Structural Characterization of Consensus:

- [TE-Aid](https://doi.org/10.1186/s13100-021-00259-7)  
  `for i in *cons.fa; do bash TE-Aid --query $i --genome <genome_assembly_name>.fna -m 400 -e 10e-10 -o TEAid_out`
- Protein domains (EMBOSS:6.6.0.0 and pfam_scan.pl tools are required)  
  `transeq -sequence <species_cdhit>.fa -outseq <species_cdhit>.fa.transeq -frame 6 -clean`  
  `pfam_scan.pl -fasta <species_cdhit>.fa.transeq -dir <pfam_database_path> -outfile <species>_pfams.txt`  

2. Phylogenetic characterization with iqtree (v 1.6.12):

Align the RVT or DDE/Tase identified within translated consensus:
`mafft --localpair --maxiterate 1000 <protein-domain.aa> > <protein-domain.aa>.maf`
Make maximum-likelihood trees:
`iqtree -s <protein-domain_seq>.maf -m MFP -nt 2 -bb 1000 -bnni`  

Sequences are renamed and classified to finally obtain the curated libraries (see article).

3. Additional filtering and classificaton.

Additional sequences might be eliminated based on the observation of the aligmnetplots (higly-fragmented/unconserved and low-copy aligments after CIAlign trimming).


2. Phylogenetic characterization.
  
A Neighbour-Joining (NJ) tree of aligned Reverse Transcriptase (RVT) domains from all the curated consensus was generated with MEGA11 (Tamura et al., 2021) to get an global overview of the types of retrotranposons.  
Based on both the NJ-tree and the structural characterization, TEs were classified as DNA transposons (class II) or as LINE, LTR or PLE (class I).





