---
Defending the PGC Specification Hypothesis
---

All of the scripts and R commands required to run our analysis and create the figures are contained within this repository. 

## Building Alignments

First off, the downloaded sequences are converted from Genbank to FASTA format if required by `extract_cds.pl`. Then 'remove_N_and_STOPs.pl` is used to remove any sequences with an ambiguous nucleotide ('N') or an internal STOP codon. These sequences are available [here](Sequences/). To identify orthologs we used BLASTX.

```
perl translate.pl Alligator_mississippiensis.cds_clean

makeblastdb -in Alligator_mississippiensis.cds_clean_prot -dbtype prot

blastx -query Alligator_sinensis.cds_clean -db Alligator_mississippiensis.cds_clean_prot -evalue 1e-06 -max_target_seqs 1 -outfmt 6 -out AS_vs_AM.txt
``` 
