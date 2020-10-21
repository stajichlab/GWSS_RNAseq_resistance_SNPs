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

#write.csv(pfam.p450s, "results/GWSS_SNP_p450s.csv")

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

