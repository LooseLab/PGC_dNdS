---
Defending the PGC Specification Hypothesis
---

All of the scripts and R commands required to run our analysis and create the figures are contained within this repository. 

## Identifying Orthologs

First off, the downloaded sequences are converted from Genbank to FASTA format if required by `extract_cds.pl`. Then `remove_N_and_STOPs.pl` is used to remove any sequences with an ambiguous nucleotide ('N') or an internal STOP codon. These sequences are available [here](Sequences/). To identify orthologs we used BLASTX, an example is shown below.

```
perl translate.pl Alligator_mississippiensis.cds_clean

makeblastdb -in Alligator_mississippiensis.cds_clean_prot -dbtype prot

blastx -query Alligator_sinensis.cds_clean -db Alligator_mississippiensis.cds_clean_prot -evalue 1e-06 -max_target_seqs 1 -outfmt 6 -out AS_vs_AM.txt
``` 
Once all pairwise BLASTs were done, we pulled out the reciprocal top hits for each pair, 3-way, 4-way and 5-way group. The orthologs that we identified are listed [here](Orthologs/).

```bash
# To pullout the pairwise and 4-way orthologs:
perl pullout_orthos.pl AM_vs_AS.txt,AS_vs_AM.txt FC_vs_FP.txt,FP_vs_FC.txt AM_vs_FP.txt,FP_vs_AM.txt

# To pullout the 3-way orthologs used for the Relative-rate test:
perl pullout_orthos_3way.pl AM_vs_FP.txt,FP_vs_AM.txt AM_vs_PT.txt,PT_vs_AM.txt 0

# To pullout the 5-way orthologs:
perl pullout_orthos_5way.pl Archosaur_4way.txt PT_vs_AM.txt,AM_vs_PT.txt
```

## Building Alignments

We first re-ran our original relative-rate tests on the independent, pairwise Falcon-Alligator-Chimp orthologs. 

```bash
# First we run the RRT using the protein sequences
perl relative_rate3_prot.pl AM_FP_PT.txt Alligator_mississippiensis.cds_clean,Falco_peregrinus.cds_clean,Pan_troglodytes.CHIMP2.1.4.cds.all.fa_clean

# We then run the relative-rate test on the whole ORF, first two codon positions and third codon position
perl relative_rate3_cod.pl AM_FP_PT.txt orf  Alligator_mississippiensis.cds_clean,Falco_peregrinus.cds_clean,Pan_troglodytes.CHIMP2.1.4.cds.all.fa_clean
perl relative_rate3_cod.pl AM_FP_PT.txt first  Alligator_mississippiensis.cds_clean,Falco_peregrinus.cds_clean,Pan_troglodytes.CHIMP2.1.4.cds.all.fa_clean
perl relative_rate3_cod.pl AM_FP_PT.txt third  Alligator_mississippiensis.cds_clean,Falco_peregrinus.cds_clean,Pan_troglodytes.CHIMP2.1.4.cds.all.fa_clean
```

Next, we recreated the data from Whittle and Extavour, running MEGA-CC on the pairs and 4-way groups of orthologs.

```bash
# We first run MEGA-CC on each pairwise and 4-way group:
perl run_megacc.pl Alligator_pairs.txt Files/ Alligator_mississippiensis.cds_clean,Alligator_sinensis.cds_clean
perl run_megacc.pl Archosaur_4way.txt Files/ Alligator_mississippiensis.cds_clean,Alligator_sinensis.cds_clean,Falco_peregrinus.cds_clean,Falco_cherrug.cds_clean

# We then pullout the dN, dS and dN/dS for the pairwise alignments:
perl get_dnds_pair.pl Alligators 15502 34469 Files/
mv summary.txt Alligator_summary.txt

# We then pullout the information from the 4-way alignments:
perl get_dnds_4way.pl Archosaurs 34470 43844 Files/
mv summary.txt Archosaur_summary.txt
```

Finally, we build the 5-way alignments and run codeml to generate the dN/dS values, and HyPhy to run the relative rate test. All of the output files from each step are included [here](Output_files/).

```bash
# We initially run the relative-rate test on the whole ORF
perl run_5way_analysis.pl 5way.txt Alligator_mississippiensis.cds_clean,Alligator_sinensis.cds_clean,Falco_peregrinus.cds_clean,Falco_cherrug.cds_clean,Pan_troglodytes.CHIMP2.1.4.cds.all.fa_clean

# We then additionally run the relative-rate test on just the first two codon positions
perl run_5way_first_rrt.pl 5way.txt Alligator_mississippiensis.cds_clean,Alligator_sinensis.cds_clean,Falco_peregrinus.cds_clean,Falco_cherrug.cds_clean,Pan_troglodytes.CHIMP2.1.4.cds.all.fa_clean
```

## Creating the figures

All of the figures were created in `R` using the following script and commands. 

```r
source("import_data.R")
library(ggplot2)
library(cowplot)
library(gtools)

individual_counts <- all_individual %>% group_by(Division) %>% summarise(lt0.5 = sum(dnds < 0.5), ut0.75 = sum(dnds  >= 0.5 & dnds < 0.75), ut1 = sum(dnds >= 0.75 & dnds < 1), gt1 = sum(dnds >= 1))
individual_counts2 <- tidyr::gather(individual_counts, "Type","Val",2:5)
individual_counts2$Type <- factor(individual_counts2$Type, levels = c("lt0.5","ut0.75","ut1","gt1"))

f1 <- ggplot(individual_counts2, aes(Division, Val, fill=Type))+geom_bar(stat="identity", position="fill")+theme(legend.position = "none")+ylab("Proportion of genes")
f2 <- ggplot(individual_counts2, aes(Division, Val, fill=Type))+geom_bar(stat="identity", position="dodge")+theme(legend.position = "none")+ylab("Number of genes")
f3 <- ggplot(all_individual, aes(Division, dnds))+geom_hline(yintercept = 0.5, linetype=2, color="lightgrey")+geom_hline(yintercept = 0.75, linetype=2, color="grey")+geom_hline(yintercept = 1, linetype=2, color="darkgrey")+geom_boxplot()
f4 <- ggplot(all_individual, aes(Division, dnds))+geom_hline(yintercept = 0.5, linetype=2, color="lightgrey")+geom_hline(yintercept = 0.75, linetype=2, color="grey")+geom_hline(yintercept = 1, linetype=2, color="darkgrey")+geom_boxplot()+coord_cartesian(ylim=c(0,1))

plot_grid(f1,f2,f3,f4)
```
