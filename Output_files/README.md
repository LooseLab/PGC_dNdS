# Output Files

The first set of files `rrtest_orf.txt`,`rrtest_first.txt`,`rrtest_third.txt` and `rrtest_prot.txt` contain information from the 3-way relative-rate test (RRT). The columns are as listed:

Column Name | Description
---|---
seq_number| arbitrary sequence identifier as in the ortholog file
prob | p-value from RRT
query_distance | distance from the common ancestor to the first species (Epigenesis)
test_distance | distance from the common ancestor to the second species (Preformation)
mam_distance | distance to the outgroup (Chimpanzee)
ln_diff | likelihood difference output from HyPhy
distinct_sites | number of parsimony-informative sites calculated in Paup
total_sites | total length of the alignment
ave_percent_identity |average percent identity (BioPerl)
overall_percent_identity | overall percent identity (BioPerl)
query_gc | %GC for the first species
test_gc | %GC for the second species
mam_gc | outgroup %GC
ave_gc | mean %GC
min_dist | minimum corrected (GTR+G) distance calculated in Paup
max_dist | maximum corrected (GTR+G) distance calculated in Paup
ave_dist | mean corrected (GTR+G) distance calculated in Paup
min_undist | minimum uncorrected distance calculated in Paup
max_undist | maxinum uncorrected distance calculated in Paup
ave_undist | mean uncorrected distance calculated in Paup

The second set of files `Alligator_summary.txt`,`Anuran_summary.txt`,`Bird_summary.txt`,`Chimp_summary.txt` and `Urodele_summary.txt` contain the pairwise dN/dS results, as well as other alignment metrics. Each file contains the following columns:

Column Name | Description
---|---
Division | group name assigned by `get_dnds_pair.pl`
SeqNumber| arbitrary sequence identifier as in the ortholog file
ID1 | first sequence ID
ID2 | second sequence ID
NoSites | total length of the alignment
dn | dN value calculated by MEGA-CC
ds | dS value calculated by MEGA-CC
dnds | dN/dS

The third set of files (`Archosaur_summary.txt`, `Amphibian_summary.txt` and `Tetrapod_summary.txt`) contain the same data but from the 4way alignments. The columns are named as such:

Column Name | Description
---|---
Division | group name assigned by `get_dnds_4way.pl`
SeqNumber| arbitrary sequence identifier as in the ortholog file
I.ID1 | first sequence ID from Epigenesis species
I.ID2 | second sequence ID from Epigenesis species
P.ID1 | first sequence ID from Preformation species
P.ID2 | second sequence ID from Preformation species
NoSites | total length of the alignment
I.dn | Epigenesis dN value calculated by MEGA-CC
I.ds | Epigenesis dS value calculated by MEGA-CC
I.dnds | Epigenesis dN/dS
P.dn | Preformation dN value calculated by MEGA-CC
P.ds | Preformation dS value calculated by MEGA-CC
P.dnds | Preformation dN/dS

The final two files contain the data from all the analyses on the 5way alignments. The first of these `5way_results.txt` contains the RRT test using all three codon positions and is as follows:

Column Name | Description
---|---
SeqNumber| arbitrary sequence identifier as in the ortholog file
id1	| species 1 (Alligator mississippiensis) ID
id2	| species 2 (Alligator sinensis) ID
id3	 | species 3 (Falco peregrinus) ID
id4	 | species 4 (Falco cherrug) ID
outgroup	| outgroup (Pan troglodytes) ID
AlnLen	| alignment length
Ave.pid	| average percent identity
Ove.pid	| overall percent identity
Ave.gc	| mean %GC
rrt_1v2.p	| RRT p-value comparing Spe1 and Spe2
rrt_1v2.d1	| RRT distance to Species 1
rrt_1v2.d2	| RRT distance to Species 2
rrt_1v3.p		| RRT p-value comparing Spe1 and Spe3
rrt_1v3.d1	| RRT distance to Species 1	
rrt_1v3.d2	| RRT distance to Species 3
rrt_1v4.p		| RRT p-value comparing Spe1 and Spe4
rrt_1v4.d1	| RRT distance to Species 1	
rrt_1v4.d2	| RRT distance to Species 4
rrt_2v3.p		| RRT p-value comparing Spe2 and Spe3
rrt_2v3.d1	| RRT distance to Species 2	
rrt_2v3.d2	| RRT distance to Species 3	
rrt_2v4.p		| RRT p-value comparing Spe2 and Spe4
rrt_2v4.d1	| RRT distance to Species 2	
rrt_2v4.d2	| RRT distance to Species 4	
rrt_3v4.p		| RRT p-value comparing Spe3 and Spe4
rrt_3v4.d1	| RRT distance to Species 3
rrt_3v4.d2		| RRT distance to Species 4
b1.dn	| AlligatorCA-Spe1 branch dN
b1.ds	| AlligatorCA-Spe1 branch dS	
b1.dnds	| AlligatorCA-Spe1 branch dN/dS	
b2.dn	| AlligatorCA-Spe2 branch dN	
b2.ds		| AlligatorCA-Spe2 branch dS
b2.dnds		| AlligatorCA-Spe2 branch dN/dS
b3.dn	| FalcoCA-Spe3 branch dN	
b3.ds	| FalcoCA-Spe3 branch dS	
b3.dnds	| FalcoCA-Spe3 branch dN/dS	
b4.dn		| FalcoCA-Spe4 branch dN
b4.ds		| FalcoCA-Spe4 branch dS
b4.dnds		| FalcoCA-Spe4 branch dN/dS
bI.dn		| Epigenesis (FalcoAlligatorCA-AlligatorCA) branch dN
bI.ds	 | Epigenesis branch dS
bI.dnds	| Epigenesis branch dN/dS
bP.dn	| Preformation branch (FalcoAlligatorCA-FalcoCA) dN
bP.ds	| Preformation branch dS
bP.dnds	| Preformation branch dN/dS
bO.dn	| Outgroup branch dN
bO.ds	| Outgroup branch dS
bO.dnds	| Outgroup branch dN/dS
p1v2.dn	| pairwise dN between Spe1 and Spe2
p1v2.ds		| pairwise dS between Spe1 and Spe2
p1v2.dnds		| pairwise dN/dS between Spe1 and Spe2
p1v3.dn		| pairwise dN between Spe1 and Spe3
p1v3.ds		| pairwise dS between Spe1 and Spe3
p1v3.dnds		| pairwise dN/dS between Spe1 and Spe3
p1v4.dn		| pairwise dN between Spe1 and Spe4
p1v4.ds		| pairwise dS between Spe1 and Spe4
p1v4.dnds	| pairwise dN/dS between Spe1 and Spe4
p2v3.dn		| pairwise dN between Spe2 and Spe3
p2v3.ds		| pairwise dS between Spe2 and Spe3
p2v3.dnds		| pairwise dN/dS between Spe2 and Spe3
p2v4.dn		| pairwise dN between Spe2 and Spe4
p2v4.ds		| pairwise dS between Spe2 and Spe4
p2v4.dnds		| pairwise dN/dS between Spe2 and Spe4
p3v4.dn		| pairwise dN between Spe3 and Spe4
p3v4.ds		| pairwise dS between Spe3 and Spe4
p3v4.dnds	| pairwise dN/dS between Spe3 and Spe4

The second 5way results file `5way_results_RRT.txt` was created using the exact same sequences as above but the relative rate test was run on only the first two codon positions. The columns are:

Column Name | Description
---|---
SeqNumber| arbitrary sequence identifier as in the ortholog file
rrtf_1v2.p	| RRT p-value comparing Spe1 and Spe2
rrtf_1v2.d1	| RRT distance to Species 1
rrtf_1v2.d2	| RRT distance to Species 2
rrtf_1v3.p		| RRT p-value comparing Spe1 and Spe3
rrtf_1v3.d1	| RRT distance to Species 1	
rrtf_1v3.d2	| RRT distance to Species 3
rrtf_1v4.p		| RRT p-value comparing Spe1 and Spe4
rrtf_1v4.d1	| RRT distance to Species 1	
rrtf_1v4.d2	| RRT distance to Species 4
rrtf_2v3.p		| RRT p-value comparing Spe2 and Spe3
rrtf_2v3.d1	| RRT distance to Species 2	
rrtf_2v3.d2	| RRT distance to Species 3	
rrtf_2v4.p		| RRT p-value comparing Spe2 and Spe4
rrtf_2v4.d1	| RRT distance to Species 2	
rrtf_2v4.d2	| RRT distance to Species 4	
rrtf_3v4.p		| RRT p-value comparing Spe3 and Spe4
rrtf_3v4.d1	| RRT distance to Species 3
rrtf_3v4.d2		| RRT distance to Species 4
