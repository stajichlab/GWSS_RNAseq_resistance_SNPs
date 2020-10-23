#GWSS SNP analysis between resistant & susceptible populations
#Pulling out all SNPs matching to p450 PFAM family
#Following GWSS_SNP_annotations.R script
#Author: Cassie Ettinger
#email: cassandra.ettinger@ucr.edu

#load libraries
library(tidyverse)
library(magrittr)

#p450 pfam family = PF00067


gwss_all_snps_annot <- read.csv("results/GWSS_SNP_annotations.snpEff.high.freqdiff.csv")

pfam.p450s <- filter(gwss_all_snps_annot, gwss_all_snps_annot$PFAM == "PF00067")

#write.csv(pfam.p450s, "results/p450/GWSS_SNP_p450s.csv")

summary(as.factor(pfam.p450s$CHROM))
#On two scaffolds
#scaffold_1338, scaffold_3544 

summary(as.factor(pfam.p450s$GENE))
#two genes
#HOVITM_016113, HOVITM_068229 

#In terminal
#head -n 1 Homalodisca_vitripennis_A6A7A9_masurca_v1.annotations.txt > GWSS_p450s_seq.txt
#grep "HOVITM_016113" HV_A6A7A9_masurca_v1.annotations.high.short.sort.txt >> GWSS_p450s_seq.txt
#grep "HOVITM_068229" HV_A6A7A9_masurca_v1.annotations.high.short.sort.txt >> GWSS_p450s_seq.txt

#opened in Excel and saved as 'csv'

#read file with seqs
p450.seqs <- read.csv("results/GWSS_p450s_seq.csv", header= TRUE)


#Make succeptible & resitant alleles 

#3 SNPs highly (diff greater than 0.5) associated with HOVITM_068229 in A vs. C comparison
#Many other SNPs with HOVITM_016113 too
#here pop A = suc & C = res        
#subset to only SNPs that differ between A/C

pfam.p450s.AC <- filter(pfam.p450s, pfam.p450s$FreqDiff_Tulare_Susceptible_A_vs_GeneralBeale_Resistant_C > 0)

#get SNP allele freqs for AvC comparison
AC.all.freq <- read.csv("intermed_results/GWSS_RNASeq1.snpEff.matrix_high.AandC_R_AlleleFreqs.csv")
AC.all.freq.clean <- subset(AC.all.freq, select = -c(A.REF.sum, C.REF.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop, X, A.vs.C.REF.freq))

#join to p450 SNP list
pfam.p450s.AC.freq <- left_join(pfam.p450s.AC, AC.all.freq.clean)

#clean up empty columns
pfam.p450s.AC.freq.clean <- subset(pfam.p450s.AC.freq, select = -c(BUSCO, Protease, CAZyme))

#save file
#write.csv(pfam.p450s.AC.freq.clean, "results/p450/GWSS_SNP_p450s.Tulare_suc_vs_Gen_B_res.csv")

#filter snps by gene
pfam.p450s.AC.HOVITM_068229 <- filter(pfam.p450s.AC.freq.clean, pfam.p450s.AC.freq.clean$GENE == "HOVITM_068229")
pfam.p450s.AC.HOVITM_016113 <- filter(pfam.p450s.AC.freq.clean, pfam.p450s.AC.freq.clean$GENE == "HOVITM_016113")


#Generate Alleles 
#Only the 3 SNPs mentioned above are different between suc & res A/C populations
pfam.p450s.AC.HOVITM_068229.v2 <- filter(pfam.p450s.AC.HOVITM_068229, pfam.p450s.AC.HOVITM_068229$FreqDiff_Susceptible_AB_vs_Resistant_CD > 0)

pfam.p450s.AC.HOVITM_068229.v2 <- pfam.p450s.AC.HOVITM_068229.v2 %>%
  mutate(SNP.pos =  readr::parse_number(CHANGEDNA))

p450.seqs.HOVITM_068229 <- filter(p450.seqs, p450.seqs$GeneID == "HOVITM_068229")

#only 3 snps so just did manually to get quick look at this gene 
#will come back to this comutationally later
#jason thought he had a script to do this 

#not doing anything with HOVITM_016113 here for now b/c no strong freq differences



#4 SNPs highly (diff greater than 0.5) associated with HOVITM_016113 in B vs. D comparison
#here pop B = suc & D = res
#subset to only SNPs that differ between A/C
#No SNPs in HOVITM_068229

pfam.p450s.BD <- filter(pfam.p450s, pfam.p450s$FreqDiff_Temecula_Susceptible_B_vs_Tulare_Resistant_D > 0)

#get SNP allele freqs for BvD comparison
BD.all.freq <- read.csv("intermed_results/GWSS_RNASeq1.snpEff.matrix_high.BandD_R_AlleleFreqs.csv")
BD.all.freq.clean <- subset(BD.all.freq, select = -c(B.REF.sum, D.REF.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop, X, B.vs.D.REF.freq))

#join to p450 SNP list
pfam.p450s.BD.freq <- left_join(pfam.p450s.BD, BD.all.freq.clean)

#clean up empty columns
pfam.p450s.BD.freq.clean <- subset(pfam.p450s.BD.freq, select = -c(BUSCO, Protease, CAZyme))

#save file
#write.csv(pfam.p450s.BD.freq.clean, "results/p450/GWSS_SNP_p450s.Temecula_suc_vs_Tulare_res.csv")
