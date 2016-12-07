# read in all data files and merge tables accordingly

library(dplyr)
library(gtools)
library(tidyr)
library(ggplot2)
library(cowplot)
library(grid)
library(gtable)

## First for the Archosaurs, we read in the Alligator, Bird and 4way datafiles. Rename some columns, and merge the data so that the 4way and pairwise alignments can be compared

arc <- read.delim("Output_files/Archosaur_summary.txt")
arc$Division = "Falco-Alligator"
bird <- read.delim("Output_files/Bird_summary.txt")
bird$Division = "Falco"
ali <- read.delim("Output_files/Alligator_summary.txt")
arc_individual <- rbind(ali, bird)

names(ali) <- c("Division", "A.SeqNumber","A.ID1","A.ID2","A.NoSites","A.dn","A.ds","A.dnds")
names(bird) <- c("Division", "B.SeqNumber","B.ID1","B.ID2","B.NoSites","B.dn","B.ds","B.dnds")

tmp1 <- inner_join(arc, ali, by = c("I.ID1" = "A.ID1"))
tmp2 <- inner_join(arc, ali, by = c("I.ID1" = "A.ID2"))
tmp2 <- plyr::rename(tmp2, replace = c("A.ID1" = "A.ID2"))
tmp <- rbind(tmp1, tmp2)
tmp1 <- inner_join(tmp, bird, by = c("P.ID1" = "B.ID1"))
tmp2 <- inner_join(tmp, bird, by = c("P.ID1" = "B.ID2"))
tmp2 <- plyr::rename(tmp2, replace = c("B.ID1" = "B.ID2"))
arc_all <- rbind(tmp1, tmp2)

## Then do the same for Amphibians

amp <- read.delim("Output_files/Amphibian_summary.txt")
amp$Division = "Xenopus'-Urodele"
anu <- read.delim("Output_files/Anuran_summary.txt")
anu$Division = "Xenopus"
uro <- read.delim("Output_files/Urodele_summary.txt")
amp_individual <- rbind(anu, uro)

names(anu) <- c("Division", "A.SeqNumber","A.ID1","A.ID2","A.NoSites","A.dn","A.ds","A.dnds")
names(uro) <- c("Division", "U.SeqNumber","U.ID1","U.ID2","U.NoSites","U.dn","U.ds","U.dnds")

tmp1 <- inner_join(amp, anu, by = c("P.ID1" = "A.ID1"))
tmp2 <- inner_join(amp, anu, by = c("P.ID1" = "A.ID2"))
tmp2 <- plyr::rename(tmp2, replace = c("A.ID1" = "A.ID2"))
tmp <- rbind(tmp1, tmp2)
tmp1 <- inner_join(tmp, uro, by = c("I.ID1" = "U.ID1"))
tmp2 <- inner_join(tmp, uro, by = c("I.ID1" = "U.ID2"))
tmp2 <- plyr::rename(tmp2, replace = c("U.ID1" = "U.ID2"))
amp_all <- rbind(tmp1, tmp2)
rm(tmp, tmp1, tmp2)

## Then for Xen vs Pan

tet <- read.delim("Output_files/Tetrapod_summary.txt")
tet$Division = "Xenopus-Pan"
anu <- read.delim("Output_files/Anuran_summary.txt")
anu$Division = "Xenopus"
pan <- read.delim("Output_files/Chimp_summary.txt")
pan$Division = "Pan"
tet_individual <- rbind(anu, pan)

names(anu) <- c("Division", "A.SeqNumber","A.ID1","A.ID2","A.NoSites","A.dn","A.ds","A.dnds")
names(pan) <- c("Division", "C.SeqNumber","C.ID1","C.ID2","C.NoSites","C.dn","C.ds","C.dnds")

tmp1 <- inner_join(tet, anu, by = c("P.ID1" = "A.ID1"))
tmp2 <- inner_join(tet, anu, by = c("P.ID1" = "A.ID2"))
tmp2 <- plyr::rename(tmp2, replace = c("A.ID1" = "A.ID2"))
tmp <- rbind(tmp1, tmp2)
tmp1 <- inner_join(tmp, pan, by = c("I.ID1" = "C.ID1"))
tmp2 <- inner_join(tmp, pan, by = c("I.ID1" = "C.ID2"))
tmp2 <- plyr::rename(tmp2, replace = c("C.ID1" = "C.ID2"))
tet_all <- rbind(tmp1, tmp2)
rm(tmp, tmp1, tmp2)

## Add extra information to each table, fold change between dN/dS values and absolute differences

arc_all$FD <- foldchange(arc_all$B.dnds, arc_all$A.dnds)
arc_all[is.na(arc_all$FD),]$FD = 0
arc_all$FDA <- foldchange(arc_all$I.dnds, arc_all$A.dnds)
arc_all[is.na(arc_all$FDA),]$FDA = 0
arc_all$FDB <- foldchange(arc_all$P.dnds, arc_all$B.dnds)
arc_all[is.na(arc_all$FDB),]$FDB = 0

amp_all$FD <- foldchange(amp_all$A.dnds, amp_all$U.dnds)
amp_all[is.na(amp_all$FD),]$FD = 0
amp_all$FDU <- foldchange(amp_all$I.dnds, amp_all$U.dnds)
amp_all[is.na(amp_all$FDU),]$FDU = 0
amp_all$FDA <- foldchange(amp_all$P.dnds, amp_all$A.dnds)
amp_all[is.na(amp_all$FDA),]$FDA = 0

tet_all$FD <- foldchange(tet_all$A.dnds, tet_all$C.dnds)
tet_all[is.na(tet_all$FD),]$FD = 0
tet_all$FDC <- foldchange(tet_all$I.dnds, tet_all$C.dnds)
tet_all[is.na(tet_all$FDC),]$FDC = 0
tet_all$FDA <- foldchange(tet_all$P.dnds, tet_all$A.dnds)
tet_all[is.na(tet_all$FDA),]$FDA = 0

arc$FD <- foldchange(arc$P.dnds, arc$I.dnds)
arc[is.na(arc$FD),]$FD = 0
tet$FD <- foldchange(tet$P.dnds, tet$I.dnds)
tet[is.na(tet$FD),]$FD = 0
amp$FD <- foldchange(amp$P.dnds, amp$I.dnds)
amp[is.na(amp$FD),]$FD = 0

arc_all$Diff <- arc_all$B.dnds - arc_all$A.dnds
amp_all$Diff <- amp_all$A.dnds - amp_all$U.dnds
tet_all$Diff <- tet_all$A.dnds - tet_all$C.dnds

arc$Diff <- arc$P.dnds - arc$I.dnds
tet$Diff <- tet$P.dnds - tet$I.dnds
amp$Diff <- amp$P.dnds - amp$I.dnds

# Get all 4way results into a table with all dN/dS values in one column (lose the gene specific information)

tmp <- data.frame(Division="Xenopus-Pan", Genus="Xenopus", ds=tet$P.ds, dn=tet$P.dn, dnds=tet$P.dnds, NoSites=tet$NoSites, SeqNumber=tet$SeqNumber )
tmp1 <- data.frame(Division="Xenopus-Pan", Genus="Pan", ds=tet$I.ds, dn=tet$I.dn, dnds=tet$I.dnds, NoSites=tet$NoSites, SeqNumber=tet$SeqNumber )
tmp2 <- data.frame(Division="Falco-Alligator", Genus="Falco", ds=arc$P.ds, dn=arc$P.dn, dnds=arc$P.dnds, NoSites=arc$NoSites, SeqNumber=arc$SeqNumber )
tmp3 <- data.frame(Division="Falco-Alligator", Genus="Alligator", ds=arc$I.ds, dn=arc$I.dn, dnds=arc$I.dnds, NoSites=arc$NoSites, SeqNumber=arc$SeqNumber )
all_4way <- rbind(tmp, tmp1, tmp2, tmp3)
all_4way$Genus <- factor(all_4way$Genus, levels=c("Falco","Alligator","Xenopus","Pan"))
rm(tmp, tmp1, tmp2, tmp3)

tmp <- data.frame(Division="Xenopus'-Urodele", Genus="Xenopus'", ds=amp$P.ds, dn=amp$P.dn, dnds=amp$P.dnds, NoSites=amp$NoSites, SeqNumber=amp$SeqNumber )
tmp1 <- data.frame(Division="Xenopus'-Urodele", Genus="Urodele", ds=amp$I.ds, dn=amp$I.dn, dnds=amp$I.dnds, NoSites=amp$NoSites, SeqNumber=amp$SeqNumber )
all_4way3 <- rbind(all_4way, tmp, tmp1)
all_4way3$Genus <- factor(all_4way3$Genus, levels=c("Falco","Alligator","Xenopus","Pan","Xenopus'","Urodele"))
all_4way3$Division <- factor(all_4way3$Division, levels=c("Falco-Alligator","Xenopus-Pan","Xenopus'-Urodele"))
rm(tmp, tmp1)

# Get all 4way results into a table with Preformation and Epigenesis dN/dS results seperated.

tet$P.ID2 <- as.factor(tet$P.ID2)
amp$P.ID2 <- as.factor(amp$P.ID2)
combi <- rbind(arc, tet)
combi$Division <- factor(combi$Division, levels=c("Falco-Alligator","Xenopus-Pan"))
combi3 <- rbind(arc,tet,amp)
combi3$Division <- factor(combi3$Division, levels=c("Falco-Alligator","Xenopus-Pan", "Xenopus'-Urodele"))

# Select the 'good' results, dS >= 0.01 and both dN and dS < 2.

good_combi <- combi %>% filter(I.ds >= 0.01, P.ds >= 0.01, I.ds < 2, P.ds < 2, I.dn < 2, P.dn < 2)
good_combi3 <- combi3 %>% filter(I.ds >= 0.01, P.ds >= 0.01, I.ds < 2, P.ds < 2, I.dn < 2, P.dn < 2)
good_tet <- tet %>% filter(I.ds >= 0.01, P.ds >= 0.01, I.ds < 2, P.ds < 2, I.dn < 2, P.dn < 2)
good_arc <- arc %>% filter(I.ds >= 0.01, P.ds >= 0.01, I.ds < 2, P.ds < 2, I.dn < 2, P.dn < 2)
good_amp <- amp %>% filter(I.ds >= 0.01, P.ds >= 0.01, I.ds < 2, P.ds < 2, I.dn < 2, P.dn < 2)

tmp <- data.frame(Division="Xenopus-Pan", Genus="Xenopus", ds=good_tet$P.ds, dn=good_tet$P.dn, dnds=good_tet$P.dnds, NoSites=good_tet$NoSites, SeqNumber=good_tet$SeqNumber )
tmp1 <- data.frame(Division="Xenopus-Pan", Genus="Pan", ds=good_tet$I.ds, dn=good_tet$I.dn, dnds=good_tet$I.dnds, NoSites=good_tet$NoSites, SeqNumber=good_tet$SeqNumber )
tmp2 <- data.frame(Division="Falco-Alligator", Genus="Falco", ds=good_arc$P.ds, dn=good_arc$P.dn, dnds=good_arc$P.dnds, NoSites=good_arc$NoSites, SeqNumber=good_arc$SeqNumber )
tmp3 <- data.frame(Division="Falco-Alligator", Genus="Alligator", ds=good_arc$I.ds, dn=good_arc$I.dn, dnds=good_arc$I.dnds, NoSites=good_arc$NoSites, SeqNumber=good_arc$SeqNumber )
good_all_4way <- rbind(tmp, tmp1, tmp2, tmp3)
good_all_4way$Genus <- factor(good_all_4way$Genus, levels=c("Falco","Alligator","Xenopus","Pan"))
rm(tmp, tmp1, tmp2, tmp3)

tmp <- data.frame(Division="Xenopus'-Urodele", Genus="Xenopus'", ds=good_amp$P.ds, dn=good_amp$P.dn, dnds=good_amp$P.dnds, NoSites=good_amp$NoSites, SeqNumber=good_amp$SeqNumber )
tmp1 <- data.frame(Division="Xenopus'-Urodele", Genus="Urodele", ds=good_amp$I.ds, dn=good_amp$I.dn, dnds=good_amp$I.dnds, NoSites=good_amp$NoSites, SeqNumber=good_amp$SeqNumber )
good_all_4way3 <- rbind(good_all_4way, tmp, tmp1)
good_all_4way3$Genus <- factor(good_all_4way3$Genus, levels=c("Falco","Alligator","Xenopus","Pan","Xenopus'","Urodele"))
good_all_4way3$Division <- factor(good_all_4way3$Division, levels=c("Falco-Alligator","Xenopus-Pan","Xenopus'-Urodele"))
rm(tmp, tmp1)

# Get the pairwise results into one table, and select the 'good' alignments

all_individual <- rbind(arc_individual, tet_individual, amp_individual[amp_individual$Division == "Urodele",])
all_individual$Division <- factor(all_individual$Division, levels=c("Falco", "Alligator","Xenopus","Pan","Urodele"))
good_all_individual <- all_individual %>% filter(ds >= 0.01, ds < 2, dn < 2)

# Parse the 3way relative rate test results, it's hard-coded to split the two independent analyses at a specific SeqNumber

rrtest_orf <- read.delim("Output_files/rrtest_orf.txt")
rrtest_first <- read.delim("Output_files/rrtest_first.txt")
rrtest_third <- read.delim("Output_files/rrtest_third.txt")
selectids <- rrtest_orf[rrtest_orf$total_sites >= 400 & !(rrtest_orf$max_dist >= 4 | rrtest_orf$max_dist-rrtest_orf$min_dist>=2) & rrtest_orf$overall_percent_identity<90,]$seq_number
rrtest <- rrtest_first[rrtest_first$seq_number %in% selectids & rrtest_first$seq_number %in% rrtest_third$seq_number,]
rrtest$Comparison = "AM_FP_PT"
rrtest[rrtest$seq_number > 9659, ]$Comparison = "AS_FC_PP"
rrtest$Result = "ns"
rrtest[rrtest$prob < 0.05 & rrtest$query_distance < rrtest$test_distance, ]$Result = "Pfast"
rrtest[rrtest$prob < 0.05 & rrtest$query_distance > rrtest$test_distance, ]$Result = "Efast"
rrt <- rrtest %>% group_by(Comparison) %>% summarise(total = n(), ns = length(Result[Result == "ns"]), Efast = length(Result[Result == "Efast"]), Pfast = length(Result[Result == "Pfast"]))
rrt$perEfast <- rrt$Efast/rrt$total*100
rrt$perPfast <- rrt$Pfast/rrt$total*100
rrt2 <- rrt %>% select(Comparison, perEfast, perPfast) %>% gather("Res","Per", 2:3)

# Parse the 5way alignments, add the branch results to get a value from the common ancestor to each species. Calculate the fold change difference between dN/dS values, and summarise the relative-rate test results.

all5 <- read.delim("Output_files/5way_results.txt")
all5$bI1.dn <- all5$b1.dn + all5$bI.dn
all5$bI1.ds <- all5$b1.ds + all5$bI.ds
all5$bI1.dnds <- all5$bI1.dn/all5$bI1.ds
all5$bI2.dn <- all5$b2.dn + all5$bI.dn
all5$bI2.ds <- all5$b2.ds + all5$bI.ds
all5$bI2.dnds <- all5$bI2.dn/all5$bI2.ds
all5$bP3.dn <- all5$b3.dn + all5$bP.dn
all5$bP3.ds <- all5$b3.ds + all5$bP.ds
all5$bP3.dnds <- all5$bP3.dn/all5$bP3.ds
all5$bP4.dn <- all5$b4.dn + all5$bP.dn
all5$bP4.ds <- all5$b4.ds + all5$bP.ds
all5$bP4.dnds <- all5$bP4.dn/all5$bP4.ds
all5$DiffbP4bI2 <- all5$bP4.dnds - all5$bI2.dnds
all5$DiffbP4bI1 <- all5$bP4.dnds - all5$bI1.dnds
all5$DiffbP3bI1 <- all5$bP3.dnds - all5$bI1.dnds
all5$DiffbP3bI2 <- all5$bP3.dnds - all5$bI2.dnds
all5$DiffbP4bI2.ds <- all5$bP4.ds - all5$bI2.ds
all5$DiffbP4bI1.ds <- all5$bP4.ds - all5$bI1.ds
all5$DiffbP3bI1.ds <- all5$bP3.ds - all5$bI1.ds
all5$DiffbP3bI2.ds <- all5$bP3.ds - all5$bI2.ds
all5$DiffbP4bI2.dn <- all5$bP4.dn - all5$bI2.dn
all5$DiffbP4bI1.dn <- all5$bP4.dn - all5$bI1.dn
all5$DiffbP3bI1.dn <- all5$bP3.dn - all5$bI1.dn
all5$DiffbP3bI2.dn <- all5$bP3.dn - all5$bI2.dn
all5$FDbP4bI2 <- foldchange(all5$bP4.dnds, all5$bI2.dnds)
all5$FDbP4bI1 <- foldchange(all5$bP4.dnds, all5$bI1.dnds)
all5$FDbP3bI1 <- foldchange(all5$bP3.dnds, all5$bI1.dnds)
all5$FDbP3bI2 <- foldchange(all5$bP3.dnds, all5$bI2.dnds)
all5$rrt_res <- "mixed"
all5[all5$rrt_1v3.p < 0.05 & all5$rrt_1v4.p < 0.05 & all5$rrt_2v3.p < 0.05 & all5$rrt_2v4.p < 0.05 & all5$rrt_1v3.d1 < all5$rrt_1v3.d2 & all5$rrt_1v4.d1 < all5$rrt_1v4.d2 & all5$rrt_2v3.d1 < all5$rrt_2v3.d2 & all5$rrt_2v4.d1 < all5$rrt_2v4.d2,]$rrt_res <- 'Pfast'
all5[all5$rrt_1v3.p < 0.05 & all5$rrt_1v4.p < 0.05 & all5$rrt_2v3.p < 0.05 & all5$rrt_2v4.p < 0.05 & all5$rrt_1v3.d1 > all5$rrt_1v3.d2 & all5$rrt_1v4.d1 > all5$rrt_1v4.d2 & all5$rrt_2v3.d1 > all5$rrt_2v3.d2 & all5$rrt_2v4.d1 > all5$rrt_2v4.d2,]$rrt_res <- 'Efast'
all5[all5$rrt_1v3.p >= 0.05 & all5$rrt_1v4.p >= 0.05 & all5$rrt_2v3.p >= 0.05 & all5$rrt_2v4.p >= 0.05,]$rrt_res <- 'ns'
all5$rrt_res <- factor(all5$rrt_res, levels=c("Pfast","ns","Efast","mixed"))

# Add the first two codon relative rate test results to this table

rrt3 <- read.delim("Output_files/5way_results_RRT.txt")
all5 <- merge(all5, rrt3, by="SeqNumber")
all5$rrtf_res <- "mixed"
all5[all5$rrtf_1v3.p < 0.05 & all5$rrtf_1v4.p < 0.05 & all5$rrtf_2v3.p < 0.05 & all5$rrtf_2v4.p < 0.05 & all5$rrtf_1v3.d1 < all5$rrtf_1v3.d2 & all5$rrtf_1v4.d1 < all5$rrtf_1v4.d2 & all5$rrtf_2v3.d1 < all5$rrtf_2v3.d2 & all5$rrtf_2v4.d1 < all5$rrtf_2v4.d2,]$rrtf_res <- 'Pfast'
all5[all5$rrtf_1v3.p < 0.05 & all5$rrtf_1v4.p < 0.05 & all5$rrtf_2v3.p < 0.05 & all5$rrtf_2v4.p < 0.05 & all5$rrtf_1v3.d1 > all5$rrtf_1v3.d2 & all5$rrtf_1v4.d1 > all5$rrtf_1v4.d2 & all5$rrtf_2v3.d1 > all5$rrtf_2v3.d2 & all5$rrtf_2v4.d1 > all5$rrtf_2v4.d2,]$rrtf_res <- 'Efast'
all5[all5$rrtf_1v3.p >= 0.05 & all5$rrtf_1v4.p >= 0.05 & all5$rrtf_2v3.p >= 0.05 & all5$rrtf_2v4.p >= 0.05,]$rrtf_res <- 'ns'
all5$rrtf_res <- factor(all5$rrtf_res, levels=c("Pfast","ns","Efast","mixed"))

##########################
# Now generate the plots #
##########################


pdf("Fig1B.pdf")
ggplot(rrt2, aes(Comparison, Per, fill=Res))+geom_bar(stat="identity", position="dodge")
dev.off()

pdf("Fig1C.pdf")
good_all4way3_means <- good_all_4way3 %>% group_by(Genus) %>% summarise(avednds = mean(dnds), aveds = mean(ds), avedn = mean(dn))
good_all4way3_means$Divergence <- c(3.7,53.0,59.8,2.8,59.8,156.2)
ggplot(good_all4way3_means, aes(Divergence, avednds, label=Genus))+geom_point()+geom_smooth(method="lm", se=F, linetype=2)+ geom_text(hjust = 0, nudge_x = 2)
dev.off()

pdf("Fig1D.pdf")
tmp <- all5[all5$rrtf_res != "mixed" & all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,] %>% select(SeqNumber, rrtf_res, bI1.dn, bI1.ds, bI1.dnds, bI2.dn, bI2.ds, bI2.dnds, bP3.dn, bP3.ds, bP3.dnds, bP4.dn, bP4.ds, bP4.dnds) %>% gather("Comparison","Value", 3:14) %>% separate(Comparison, c("Spe","test")) %>% spread(test, Value)
tmp[tmp$Spe == "bI1",]$Spe = "AM"
tmp[tmp$Spe == "bI2",]$Spe = "AS"
tmp[tmp$Spe == "bP3",]$Spe = "FP"
tmp[tmp$Spe == "bP4",]$Spe = "FC"
tmp$Spe = factor(tmp$Spe, levels=c("AM","AS","FP","FC"))
ggplot(tmp, aes(rrtf_res, dnds, color=Spe))+geom_boxplot(outlier.shape = NA)+coord_cartesian(ylim=c(0,0.5))
dev.off()

pdf("FigS2.pdf")
individual_counts <- all_individual %>% group_by(Division) %>% summarise(lt0.5 = sum(dnds < 0.5), ut0.75 = sum(dnds  >= 0.5 & dnds < 0.75), ut1 = sum(dnds >= 0.75 & dnds < 1), gt1 = sum(dnds >= 1))
individual_counts2 <- tidyr::gather(individual_counts, "Type","Val",2:5)
individual_counts2$Type <- factor(individual_counts2$Type, levels = c("lt0.5","ut0.75","ut1","gt1"))
f1 <- ggplot(individual_counts2, aes(Division, Val, fill=Type))+geom_bar(stat="identity", position="fill")+theme(legend.position = "none")
f2 <- ggplot(individual_counts2, aes(Division, Val, fill=Type))+geom_bar(stat="identity", position="dodge")+theme(legend.position = "none")
f3 <- ggplot(all_individual, aes(Division, dnds))+geom_hline(yintercept = 0.5, linetype=2, color="lightgrey")+geom_hline(yintercept = 0.75, linetype=2, color="grey")+geom_hline(yintercept = 1, linetype=2, color="darkgrey")+geom_boxplot()
f4 <- ggplot(all_individual, aes(Division, dnds))+geom_hline(yintercept = 0.5, linetype=2, color="lightgrey")+geom_hline(yintercept = 0.75, linetype=2, color="grey")+geom_hline(yintercept = 1, linetype=2, color="darkgrey")+geom_boxplot()+coord_cartesian(ylim=c(0,1))
plot_grid(f1,f2,f3,f4)
dev.off()

pdf("FigS3.pdf")
f1 <- ggplot(arc_all, aes(A.dnds, I.dnds))+geom_abline(intercept=0)+geom_point(size=0.5)+geom_point(data=arc_all[arc_all$FDA >= 1.5 | arc_all$FDA <= -1.5,], color="blue", size=0.5)+xlab("dN/dS pairwise alignment")+ylab("dN/dS 4way alignment")+ggtitle("Alligator")+theme(axis.text = element_text(size=8), axis.title = element_text(size=10), plot.title = element_text(size=12))
f2 <- ggplot(arc_all, aes(B.dnds, P.dnds))+geom_abline(intercept=0)+geom_point(size=0.5)+geom_point(data=arc_all[arc_all$FDB >= 1.5 | arc_all$FDB <= -1.5,], color="blue", size=0.5)+xlab("dN/dS pairwise alignment")+ylab("dN/dS 4way alignment")+ggtitle("Falco")+theme(axis.text = element_text(size=8), axis.title = element_text(size=10), plot.title = element_text(size=12))
f3 <- ggplot(tet_all, aes(A.dnds, P.dnds))+geom_abline(intercept=0)+geom_point(size=0.5)+geom_point(data=tet_all[tet_all$FDA >= 1.5 | tet_all$FDA <= -1.5,], color="blue", size=0.5)+xlab("dN/dS pairwise alignment")+ylab("dN/dS 4way alignment")+ggtitle("Xenopus")+theme(axis.text = element_text(size=8), axis.title = element_text(size=10), plot.title = element_text(size=12))
f4 <- ggplot(tet_all, aes(C.dnds, I.dnds))+geom_abline(intercept=0)+geom_point(size=0.5)+geom_point(data=tet_all[tet_all$FDC >= 1.5 | tet_all$FDC <= -1.5,], color="blue", size=0.5)+xlab("dN/dS pairwise alignment")+ylab("dN/dS 4way alignment")+ggtitle("Pan")+theme(axis.text = element_text(size=8), axis.title = element_text(size=10), plot.title = element_text(size=12))
f5 <- ggplot(amp_all, aes(A.dnds, P.dnds))+geom_abline(intercept=0)+geom_point(size=0.5)+geom_point(data=amp_all[amp_all$FDA >= 1.5 | amp_all$FDA <= -1.5,], color="blue", size=0.5)+xlab("dN/dS pairwise alignment")+ylab("dN/dS 4way alignment")+ggtitle("Xenopus'")+theme(axis.text = element_text(size=8), axis.title = element_text(size=10), plot.title = element_text(size=12))
f6 <- ggplot(amp_all, aes(U.dnds, I.dnds))+geom_abline(intercept=0)+geom_point(size=0.5)+geom_point(data=amp_all[amp_all$FDU >= 1.5 | amp_all$FDU <= -1.5,], color="blue", size=0.5)+xlab("dN/dS pairwise alignment")+ylab("dN/dS 4way alignment")+ggtitle("Urodele")+theme(axis.text = element_text(size=8), axis.title = element_text(size=10), plot.title = element_text(size=12))
plot_grid(f2,f1,f3,f4,f5,f6, nrow=3)
dev.off()

pdf("FigS4.pdf")
f1 <- ggplot(all_4way3, aes(ds, dnds, color=Genus))+geom_point(size=0.5)+coord_cartesian(xlim=c(0,2), ylim=c(0,10.8))+theme(legend.position=c(0.85,0.85), legend.title = element_blank(), plot.title = element_text(hjust=0))+ggtitle("a")+xlab("dS")+ylab("dN/dS")
f2 <- ggplot(all_4way3[all_4way3$ds >= 0.01,], aes(ds, dnds, color=Genus))+geom_point(size=0.5)+coord_cartesian(xlim=c(0,2), ylim=c(0,10.8))+theme(legend.position="none", plot.title = element_text(hjust=0))+ggtitle("b")+xlab("dS")+ylab("dN/dS")
plot_grid(f1,f2, nrow=2)
dev.off()

pdf("FigS5.pdf")
good_Fway_counts <- good_all_4way3 %>% group_by(Genus) %>% summarise(lt0.5 = sum(dnds < 0.5), ut0.75 = sum(dnds  >= 0.5 & dnds < 0.75), ut1 = sum(dnds >= 0.75 & dnds < 1), gt1 = sum(dnds >= 1))
good_Fway_counts2 <- tidyr::gather(good_Fway_counts, "Type","Val",2:5)
good_Fway_counts2$Type <- factor(good_Fway_counts2$Type, levels = c("lt0.5","ut0.75","ut1","gt1"))
f1 <- ggplot(good_Fway_counts2, aes(Genus, Val, fill=Type))+geom_bar(stat="identity", position="fill")+theme(legend.position = "none")
f2 <- ggplot(good_Fway_counts2, aes(Genus, Val, fill=Type))+geom_bar(stat="identity", position="dodge")+theme(legend.position = "none")
f3 <- ggplot(good_all_4way3, aes(Genus, dnds))+geom_hline(yintercept = 0.5, linetype=2, color="lightgrey")+geom_hline(yintercept = 0.75, linetype=2, color="grey")+geom_hline(yintercept = 1, linetype=2, color="darkgrey")+geom_boxplot()
f4 <- ggplot(good_all_4way3, aes(Genus, dnds))+geom_hline(yintercept = 0.5, linetype=2, color="lightgrey")+geom_hline(yintercept = 0.75, linetype=2, color="grey")+geom_hline(yintercept = 1, linetype=2, color="darkgrey")+geom_boxplot()+coord_cartesian(ylim=c(0,1))
plot_grid(f1,f2,f3,f4)
dev.off()

pdf("FigS6.pdf")
f1 <- ggplot(all5)+geom_freqpoly(aes(bI1.ds),binwidth=0.01, color="blue")+geom_freqpoly(aes(bI2.ds),binwidth=0.01, color="skyblue")+geom_freqpoly(aes(bP3.ds),binwidth=0.01, color="red")+geom_freqpoly(aes(bP4.ds),binwidth=0.01, color="pink")+coord_cartesian(xlim=c(0,1))+geom_vline(xintercept = 0.01, linetype=2, color="grey")+geom_vline(xintercept = 0.04, linetype=2)
f2 <- ggplot(all5)+geom_freqpoly(aes(bI1.dn),binwidth=0.001, color="blue")+geom_freqpoly(aes(bI2.dn),binwidth=0.001, color="skyblue")+geom_freqpoly(aes(bP3.dn),binwidth=0.001, color="red")+geom_freqpoly(aes(bP4.dn),binwidth=0.001, color="pink")+coord_cartesian(xlim=c(0,0.2))
plot_grid(f1,f2, nrow=2)
dev.off()

pdf("FigS7.pdf")
f1 <- ggplot(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,])+geom_boxplot(aes("AM",bI1.dnds), color="blue")+geom_boxplot(aes("AS", bI2.dnds), color="skyblue")+geom_boxplot(aes("FP", bP3.dnds), color="red")+geom_boxplot(aes("FC", bP4.dnds), color="pink")
f2 <- ggplot(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,])+geom_boxplot(aes("AM",bI1.dnds), color="blue", outlier.shape = NA)+geom_boxplot(aes("AS", bI2.dnds), color="skyblue", outlier.shape = NA)+geom_boxplot(aes("FP", bP3.dnds), color="red", outlier.shape = NA)+geom_boxplot(aes("FC", bP4.dnds), color="pink", outlier.shape = NA)+coord_cartesian(ylim=c(0,0.4))
plot_grid(f1,f2)
dev.off()

pdf("FigS8.pdf", onefile = F)
p1 <- ggplot(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,], aes(bI1.dnds, bP3.dnds))+geom_abline()+geom_point(size=1)+geom_point(data = all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2 & (all5$FDbP3bI1 >= 1.5 | all5$FDbP3bI1 <= -1.5),], color="blue", size=1)+xlab("Alligator mississippiensis dN/dS")+ylab("Falco peregrinus dN/dS")+scale_x_continuous(expand=c(0,0))+scale_y_continuous(expand=c(0,0))+expand_limits(y=c(min(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,]$bP3.dnds) - 0.1, max(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,]$bP3.dnds) + 0.1))+expand_limits(x = c(min(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,]$bI1.dnds) - 0.1, max(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,]$bI1.dnds) + 0.1))+ theme_grey()+theme(plot.margin = unit(c(0.2,0.2,0.5,0.5), "lines"))
p2 <- ggplot(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,], aes(x = factor(1), y = bI1.dnds))+geom_boxplot(outlier.shape = 1, outlier.size = 1)+scale_y_continuous(expand=c(0,0))+expand_limits(y=c(min(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,]$bI1.dnds) - 0.1, max(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,]$bI1.dnds) + 0.1))+coord_flip()+ theme_grey()+theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(1,0.2,-0.5,0.5), "lines"))
p3 <- ggplot(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,], aes(x = factor(1), y = bP3.dnds))+geom_boxplot(outlier.shape = 1, outlier.size = 1)+scale_y_continuous(expand=c(0,0))+expand_limits(y=c(min(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,]$bP3.dnds) - 0.1, max(all5[all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,]$bP3.dnds) + 0.1))+ theme_grey()+theme(axis.text = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), plot.margin = unit(c(0.2,1,0.5,-0.5), "lines"))
gt1 <- ggplot_gtable(ggplot_build(p1))
gt2 <- ggplot_gtable(ggplot_build(p2))
gt3 <- ggplot_gtable(ggplot_build(p3))
maxWidth <- unit.pmax(gt1$widths[2:3], gt2$widths[2:3])
maxHeight <- unit.pmax(gt1$heights[4:5], gt3$heights[4:5])
gt1$widths[2:3]  <- as.list(maxWidth)
gt2$widths[2:3]  <- as.list(maxWidth)
gt1$heights[4:5] <- as.list(maxHeight)
gt3$heights[4:5] <- as.list(maxHeight)
gt <- gtable(widths = unit(c(7,1), "null"), height = unit(c(1,7), "null"))
gt <- gtable_add_grob(gt, gt1, 2, 1)
gt <- gtable_add_grob(gt, gt2, 1, 1)
gt <- gtable_add_grob(gt, gt3, 2, 2)
grid.newpage()
grid.draw(gt)
dev.off()

pdf("FigS9.pdf")
tmp <- all5[all5$rrtf_res != "mixed" & all5$bI1.ds > 0.04 & all5$bI2.ds > 0.04 & all5$bP3.ds > 0.04 & all5$bP4.ds > 0.04 & all5$bP3.ds < 2 & all5$bI1.ds < 2 & all5$bP4.ds < 2 & all5$bI2.ds < 2 & all5$bP3.dn < 2 & all5$bI1.dn < 2 & all5$bP4.dn < 2 & all5$bI2.dn < 2,] %>% select(SeqNumber, rrtf_res, bI1.dn, bI1.ds, bI1.dnds, bI2.dn, bI2.ds, bI2.dnds, bP3.dn, bP3.ds, bP3.dnds, bP4.dn, bP4.ds, bP4.dnds) %>% gather("Comparison","Value", 3:14) %>% separate(Comparison, c("Spe","test")) %>% spread(test, Value)
tmp[tmp$Spe == "bI1",]$Spe = "AM"
tmp[tmp$Spe == "bI2",]$Spe = "AS"
tmp[tmp$Spe == "bP3",]$Spe = "FP"
tmp[tmp$Spe == "bP4",]$Spe = "FC"
tmp$Spe = factor(tmp$Spe, levels=c("AM","AS","FP","FC"))
f1 <- ggplot(tmp, aes(rrtf_res, dnds, color=Spe))+geom_boxplot(outlier.shape = NA)+coord_cartesian(ylim=c(0,0.5))
f2 <- ggplot(tmp, aes(rrtf_res, dn, color=Spe))+geom_boxplot(outlier.shape = NA)+coord_cartesian(ylim=c(0,0.25))
f3 <- ggplot(tmp, aes(rrtf_res, ds, color=Spe))+geom_boxplot(outlier.shape = NA)+coord_cartesian(ylim=c(0,1))
plot_grid(f1,f2,f3, nrow=3)
dev.off()

pdf("FigS10.pdf")
tmp <- all5[all5$p1v2.ds < 2 & all5$p1v2.ds > 0.01 & all5$p3v4.ds < 2 & all5$p3v4.ds > 0.01,] %>% select(SeqNumber, rrtf_res, p1v2.dn, p1v2.ds, p1v2.dnds, p3v4.ds, p3v4.dn, p3v4.dnds) %>% gather("Comparison","Value", 3:8) %>% separate(Comparison, c("Spe","test")) %>% spread(test, Value)
tmp[tmp$Spe == "p1v2",]$Spe = "Alligator"
tmp[tmp$Spe == "p3v4",]$Spe = "Falco"
f1 <- ggplot(tmp[tmp$rrtf_res != "mixed",], aes(rrtf_res, dnds, color=Spe))+geom_boxplot(outlier.shape = NA)+coord_cartesian(ylim=c(0,1.5))
f2 <- ggplot(tmp[tmp$rrtf_res != "mixed",], aes(rrtf_res, dn, color=Spe))+geom_boxplot(outlier.shape = NA)+coord_cartesian(ylim=c(0,0.03))
f3 <- ggplot(tmp[tmp$rrtf_res != "mixed",], aes(rrtf_res, ds, color=Spe))+geom_boxplot(outlier.shape = NA)+coord_cartesian(ylim=c(0,0.125))
plot_grid(f1,f2,f3, nrow=3)
dev.off()