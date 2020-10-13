#GWSS SNP analysis between resistant & susceptible populations
#Only comparing A vs. C
#Author: Cassie Ettinger
#email: cassandra.ettinger@ucr.edu

#load libraries
library(tidyverse)
library(magrittr)

#Remove B & D columns from files to save space & save in terminal 
#also! R was not happy about empty columns (??) so this solved that issue
# cat GWSS_RNASeq1.snpEff.matrix.tab | cut -f 1-14,19-22,27 > GWSS_RNASeq1.snpEff.matrix.AandC.tab 
# cat GWSS_RNASeq1.snpEff.matrix_high.tab | cut -f 1-14,19-22,27 > GWSS_RNASeq1.snpEff.matrix_high.AandC.tab 


#load in snp data from Jason (this is after cutting out B/D pops)
#note here we are starting with 'high' table which Jason said was filtered in some way
snps <- read.delim("data/GWSS_RNASeq1.snpEff.matrix_high.AandC.tab")


#remove brackets around alternate alleles
for (i in 1:length(snps$ALT)) {
  snps$ALT[i] <- str_sub(snps$ALT[i], 2, -2)
}


#make list of all our individuals
pop_col = c('A.1', 'A.2', 'A.3', 'A.4', 'C.1', 'C.2', 'C.3', 'C.4')


#remove trailing '/' from populations
for (pop in pop_col) {
  for (i in 1:length(snps$GENE)) {
    snps[,pop][i] <- str_sub(snps[,pop][i], 0, -2)
  }
}


#subset to only be SNPs (no insertions / deletions) 
#note the way I have done this means that if there are 2 alternate alleles (E.g. A, G)
#those would not be included here - but most of those scenarios involve the "*" (deletion) relative to ref so prob OK

#double checking what I just said to make sure 
summary(as.factor(snps$ALT))


#add column to indicate if row contains indel
snps %<>%
  mutate(indel = ifelse(nchar(ALT) > 1 | nchar(REF) > 1, TRUE, FALSE))


#filter to only be SNPs, can change indel -> FALSE to get insertions / deletions relative to REF instead         
snps.filt <- filter(snps, indel == FALSE)


#filter to only get SNPs sequenced across all individuals 
#ask Jason if we want to pursue SNPs where some individuals were not covered

#set up column with all data assumed to be present
snps.filt$missing <- FALSE

#for each individual, check if SNP called, if missing - reassign value to TRUE
#else leave value as is, with default of FALSE
for (pop in pop_col) {
  
  snps.filt %<>%
    mutate(missing = ifelse(snps.filt[pop] == ".", TRUE, as.character(missing)))
  
}

#filter to only have table where SNPS were covered for all individuals
snps.filt.cov <- filter(snps.filt, missing == FALSE)


# /:  genotype unphased◦ |:  genotype phased
# I think we can ignore these for now as we only want allele freqs

#calculate allele frequencies for reference allele 
#b/c of the way we filtered - I think we should only see the REF & ALT nucleotides

#set up columns to save sums
snps.filt.cov$A.REF.sum <- 0
snps.filt.cov$C.REF.sum <- 0

#set up indivduals per popualtion
A.pop <- c("A.1", "A.2", "A.3", "A.4")
C.pop <- c("C.1", "C.2", "C.3", "C.4")

#sum A pop ref alleles
for (pop in A.pop) {
  for (i in 1:length(snps.filt.cov$GENE)) {
    snps.filt.cov$A.REF.sum[i] <- snps.filt.cov$A.REF.sum[i] + str_count(snps.filt.cov[,pop][i], pattern = snps.filt.cov$REF[i])
  }
}

#sum C pop ref alleles
for (pop in C.pop) {
  for (i in 1:length(snps.filt.cov$GENE)) {
    snps.filt.cov$C.REF.sum[i] <- snps.filt.cov$C.REF.sum[i] + str_count(snps.filt.cov[,pop][i], pattern = snps.filt.cov$REF[i])
  }
}


#set up to save freq
snps.filt.cov$A.REF.freq <- 0
snps.filt.cov$C.REF.freq <- 0
snps.filt.cov$A.ALT.freq <- 0
snps.filt.cov$C.ALT.freq <- 0

#turn sums to proportion & calculate alt allele freq based on REF
for (i in 1:length(snps.filt.cov$GENE)) {
  snps.filt.cov$A.REF.freq[i] <- (snps.filt.cov$A.REF.sum[i] / (2*length(A.pop))) #2 times length bc diploid
  snps.filt.cov$C.REF.freq[i] <- (snps.filt.cov$C.REF.sum[i] / (2*length(C.pop))) #2 times length bc diploid
  
  snps.filt.cov$A.ALT.freq[i] <- 1 - snps.filt.cov$A.REF.freq[i]
  snps.filt.cov$C.ALT.freq[i] <- 1 - snps.filt.cov$C.REF.freq[i]
}


snps.filt.cov$A.vs.C.REF.freq <- 0

#calculate abs difference in allele freqs between populations
for (i in 1:length(snps.filt.cov$GENE)) {
  snps.filt.cov$A.vs.C.REF.freq[i] <- abs(as.numeric(snps.filt.cov$A.REF.freq[i]) - as.numeric(snps.filt.cov$C.REF.freq[i])) 
}


#filter out alleles where there is no difference
snps.filt.cov.var <- filter(snps.filt.cov, A.vs.C.REF.freq != 0)


#order by largest to smallest difference
snps.filt.cov.var.sort <- snps.filt.cov.var[order(-snps.filt.cov.var$A.vs.C.REF.freq),]

#Visualize difference in freqs
hist(snps.filt.cov.var.sort$A.vs.C.REF.freq, xlab = "Absolute difference between SNP frequencies for reference allele", main = "Population A (Tulare, suc) vs. Population C (General Beale, res)")

#visualize by base
ggplot(snps.filt.cov.var.sort, aes(x=A.vs.C.REF.freq)) +geom_histogram() + facet_grid(~REF) +xlab("Absolute difference between SNP frequencies for reference allele") + ggtitle("Population A (Tulare, suc) vs. Population C (General Beale, res)")
#ggsave(filename = 'plots/A.Tulare.Suc.vs.C.GeneralBeale.Res.Ref.Freq.pdf', plot = last_plot(), device = 'pdf', width = 8, height = 4, dpi = 300)


#Visualize individual freqs & by base
ggplot(snps.filt.cov.var.sort, aes(x=A.REF.freq)) +geom_histogram() + facet_grid(~REF) +xlab("SNP frequences for reference allele") + ggtitle("Population A (Tulare, suc)")
#ggsave(filename = 'plots/A.Tulare.Suc.Ref.Freq.pdf', plot = last_plot(), device = 'pdf', width = 8, height = 4, dpi = 300)


ggplot(snps.filt.cov.var.sort, aes(x=C.REF.freq)) +geom_histogram() + facet_grid(~REF) +xlab("SNP frequences for reference allele") + ggtitle("Population C (General Beale, res)")
#ggsave(filename = 'plots/C.GeneralBeale.Res.Ref.Freq.pdf', plot = last_plot(), device = 'pdf', width = 8, height = 4, dpi = 300)


#save file as csv
#write.csv(snps.filt.cov.var.sort, file = "results/GWSS_RNASeq1.snpEff.matrix_high.AandC_R_AlleleFreqs.csv")

#make short file
#remove some columns generated in R
snps.filt.cov.var.sort.short <- subset(snps.filt.cov.var.sort, select = -c(A.REF.sum, C.REF.sum, indel, missing))

#filter out alleles that don't differ > .5
snps.filt.cov.var.sort.short <- filter(snps.filt.cov.var.sort.short, A.vs.C.REF.freq > 0.5)

#save file as csv
#write.csv(snps.filt.cov.var.sort.short, file = "results/GWSS_RNASeq1.snpEff.matrix_high.AandC_R_AlleleFreqs.short.csv")




