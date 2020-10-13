#GWSS SNP analysis between resistant & susceptible populations
#By location & resitance status
#Author: Cassie Ettinger
#email: cassandra.ettinger@ucr.edu

#load libraries
library(tidyverse)
library(magrittr)


#load in snp data from Jason 
#note here we are starting with 'high' table which Jason said was filtered in some way
#used cut to only select columns 1-27 (extra trailing empty columns angered R gods)
snps <- read.delim("data/GWSS_RNASeq1.snpEff.matrix_high.fix.tab")


#remove brackets around alternate alleles
for (i in 1:length(snps$ALT)) {
  snps$ALT[i] <- str_sub(snps$ALT[i], 2, -2)
}


#make list of all our individuals
pop_col = c('A.1', 'A.2', 'A.3', 'A.4', 
            'B.1', 'B.2', 'B.3', 'B.4', 
            'C.1', 'C.2', 'C.3', 'C.4', 
            'D.1', 'D.2','D.3','D.4')


#set up indivduals per popualtion
A.pop <- c("A.1", "A.2", "A.3", "A.4")
B.pop <- c("B.1", "B.2", "B.3", "B.4")
C.pop <- c("C.1", "C.2", "C.3", "C.4")
D.pop <- c("D.1", "D.2", "D.3", "D.4")


#remove trailing '/' from populations
#this takes time & prob is not necessary
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
snps.filt$missing.Apop <- FALSE
snps.filt$missing.Bpop <- FALSE
snps.filt$missing.Cpop <- FALSE
snps.filt$missing.Dpop <- FALSE

#for each individual, check if SNP called, if missing - reassign value to TRUE
#else leave value as is, with default of FALSE
for (pop in pop_col) {
  
  snps.filt %<>%
    mutate(missing = ifelse(snps.filt[pop] == ".", TRUE, as.character(missing)))
  
  if (pop %in% A.pop) {
    snps.filt %<>%
      mutate(missing.Apop = ifelse(snps.filt[pop] == ".", TRUE, as.character(missing.Apop)))
  } else if (pop %in% B.pop) {
      snps.filt %<>%
        mutate(missing.Bpop = ifelse(snps.filt[pop] == ".", TRUE, as.character(missing.Bpop)))
    
       } else if (pop %in% C.pop) {
           snps.filt %<>%
        mutate(missing.Cpop = ifelse(snps.filt[pop] == ".", TRUE, as.character(missing.Cpop)))
      
           } else {
              snps.filt %<>%
          mutate(missing.Dpop = ifelse(snps.filt[pop] == ".", TRUE, as.character(missing.Dpop)))
        
             }
}

#filter to only have table where SNPS were covered for all individuals across all pops
snps.filt.cov <- filter(snps.filt, missing == FALSE)

#filter to only have table where SNPS were covered for all individuals in comparison being made
snps.filt.cov.AC <- filter(snps.filt, missing.Apop == FALSE &  missing.Cpop == FALSE)
snps.filt.cov.AD <- filter(snps.filt, missing.Apop == FALSE &  missing.Dpop == FALSE)
snps.filt.cov.BC <- filter(snps.filt, missing.Bpop == FALSE &  missing.Cpop == FALSE)
snps.filt.cov.BD <- filter(snps.filt, missing.Bpop == FALSE &  missing.Dpop == FALSE)

# /:  genotype unphased◦ |:  genotype phased
# I think we can ignore these for now as we only want allele freqs

#calculate allele frequencies for reference allele 
#b/c of the way we filtered - I think we should only see the REF & ALT nucleotides

#set up columns to save sums for whole dataset
snps.filt.cov$A.REF.sum <- 0
snps.filt.cov$B.REF.sum <- 0
snps.filt.cov$C.REF.sum <- 0
snps.filt.cov$D.REF.sum <- 0

#set up columns to save sums for each pop comparison
snps.filt.cov.AC$A.REF.sum <- 0
snps.filt.cov.AC$C.REF.sum <- 0

snps.filt.cov.AD$A.REF.sum <- 0
snps.filt.cov.AD$D.REF.sum <- 0

snps.filt.cov.BC$B.REF.sum <- 0
snps.filt.cov.BC$C.REF.sum <- 0

snps.filt.cov.BD$B.REF.sum <- 0
snps.filt.cov.BD$D.REF.sum <- 0



#sum A pop ref alleles
for (pop in A.pop) {
  for (i in 1:length(snps.filt.cov$GENE)) {
    snps.filt.cov$A.REF.sum[i] <- snps.filt.cov$A.REF.sum[i] + str_count(snps.filt.cov[,pop][i], pattern = snps.filt.cov$REF[i])
  }
  for (i in 1:length(snps.filt.cov.AC$GENE)) {
    snps.filt.cov.AC$A.REF.sum[i] <- snps.filt.cov.AC$A.REF.sum[i] + str_count(snps.filt.cov.AC[,pop][i], pattern = snps.filt.cov.AC$REF[i])
  }
  for (i in 1:length(snps.filt.cov.AD$GENE)) {
    snps.filt.cov.AD$A.REF.sum[i] <- snps.filt.cov.AD$A.REF.sum[i] + str_count(snps.filt.cov.AD[,pop][i], pattern = snps.filt.cov.AD$REF[i])
  }
}

#sum B pop ref alleles
for (pop in B.pop) {
  for (i in 1:length(snps.filt.cov$GENE)) {
    snps.filt.cov$B.REF.sum[i] <- snps.filt.cov$B.REF.sum[i] + str_count(snps.filt.cov[,pop][i], pattern = snps.filt.cov$REF[i])
  }
  for (i in 1:length(snps.filt.cov.BC$GENE)) {
    snps.filt.cov.BC$B.REF.sum[i] <- snps.filt.cov.BC$B.REF.sum[i] + str_count(snps.filt.cov.BC[,pop][i], pattern = snps.filt.cov.BC$REF[i])
  }
  for (i in 1:length(snps.filt.cov.BD$GENE)) {
    snps.filt.cov.BD$B.REF.sum[i] <- snps.filt.cov.BD$B.REF.sum[i] + str_count(snps.filt.cov.BD[,pop][i], pattern = snps.filt.cov.BD$REF[i])
  }
}

#sum C pop ref alleles
for (pop in C.pop) {
  for (i in 1:length(snps.filt.cov$GENE)) {
    snps.filt.cov$C.REF.sum[i] <- snps.filt.cov$C.REF.sum[i] + str_count(snps.filt.cov[,pop][i], pattern = snps.filt.cov$REF[i])
  }
  for (i in 1:length(snps.filt.cov.BC$GENE)) {
    snps.filt.cov.BC$C.REF.sum[i] <- snps.filt.cov.BC$C.REF.sum[i] + str_count(snps.filt.cov.BC[,pop][i], pattern = snps.filt.cov.BC$REF[i])
  }
  for (i in 1:length(snps.filt.cov.AC$GENE)) {
    snps.filt.cov.AC$C.REF.sum[i] <- snps.filt.cov.AC$C.REF.sum[i] + str_count(snps.filt.cov.AC[,pop][i], pattern = snps.filt.cov.AC$REF[i])
  }
}


#sum D pop ref alleles
for (pop in D.pop) {
  for (i in 1:length(snps.filt.cov$GENE)) {
    snps.filt.cov$D.REF.sum[i] <- snps.filt.cov$D.REF.sum[i] + str_count(snps.filt.cov[,pop][i], pattern = snps.filt.cov$REF[i])
  }
  for (i in 1:length(snps.filt.cov.BD$GENE)) {
    snps.filt.cov.BD$D.REF.sum[i] <- snps.filt.cov.BD$D.REF.sum[i] + str_count(snps.filt.cov.BD[,pop][i], pattern = snps.filt.cov.BD$REF[i])
  }
  for (i in 1:length(snps.filt.cov.AD$GENE)) {
    snps.filt.cov.AD$D.REF.sum[i] <- snps.filt.cov.AD$D.REF.sum[i] + str_count(snps.filt.cov.AD[,pop][i], pattern = snps.filt.cov.AD$REF[i])
  }
}

#A - Tulare (suc)
#B - Temecula (suc)
#C - GeneralBeale (res)
#D - Tulare (res)


#set up to save lumped by resistance status
snps.filt.cov$A.B.suc.REF.freq <- 0
snps.filt.cov$C.D.res.REF.freq <-0

snps.filt.cov$A.B.suc.ALT.freq <- 0
snps.filt.cov$C.D.res.ALT.freq <-0


#set up to save freq for each pop 
snps.filt.cov.AC$A.REF.freq <- 0
snps.filt.cov.AC$C.REF.freq <- 0

snps.filt.cov.AD$A.REF.freq <- 0
snps.filt.cov.AD$D.REF.freq <- 0

snps.filt.cov.BC$B.REF.freq <- 0
snps.filt.cov.BC$C.REF.freq <- 0

snps.filt.cov.BD$B.REF.freq <- 0
snps.filt.cov.BD$D.REF.freq <- 0

#alt freqs
snps.filt.cov.AC$A.ALT.freq <- 0
snps.filt.cov.AC$C.ALT.freq <- 0

snps.filt.cov.AD$A.ALT.freq <- 0
snps.filt.cov.AD$D.ALT.freq <- 0

snps.filt.cov.BC$B.ALT.freq <- 0
snps.filt.cov.BC$C.ALT.freq <- 0

snps.filt.cov.BD$B.ALT.freq <- 0
snps.filt.cov.BD$D.ALT.freq <- 0


#turn sums to proportion & calculate alt allele freq based on REF
for (i in 1:length(snps.filt.cov$GENE)) {

  snps.filt.cov$A.B.suc.REF.freq[i] <- ((snps.filt.cov$A.REF.sum[i]+snps.filt.cov$B.REF.sum[i]) / (2*(length(A.pop)+length(B.pop)))) 
  snps.filt.cov$C.D.res.REF.freq[i] <- ((snps.filt.cov$C.REF.sum[i]+snps.filt.cov$D.REF.sum[i]) / (2*(length(C.pop)+length(D.pop)))) 
  
  snps.filt.cov.AC$A.REF.freq[i] <- (snps.filt.cov.AC$A.REF.sum[i] / (2*length(A.pop))) #2 times length bc diploid
  snps.filt.cov.AC$C.REF.freq [i] <- (snps.filt.cov.AC$C.REF.sum[i] / (2*length(C.pop))) #2 times length bc diploid

  snps.filt.cov.AD$A.REF.freq[i] <- (snps.filt.cov.AD$A.REF.sum[i] / (2*length(A.pop))) #2 times length bc diploid
  snps.filt.cov.AD$D.REF.freq [i] <- (snps.filt.cov.AD$D.REF.sum[i] / (2*length(D.pop))) #2 times length bc diploid
  
  snps.filt.cov.BC$B.REF.freq[i] <- (snps.filt.cov.BC$B.REF.sum[i] / (2*length(B.pop))) #2 times length bc diploid
  snps.filt.cov.BC$C.REF.freq [i] <- (snps.filt.cov.BC$C.REF.sum[i] / (2*length(C.pop))) #2 times length bc diploid
  
  snps.filt.cov.BD$B.REF.freq[i] <- (snps.filt.cov.BD$B.REF.sum[i] / (2*length(B.pop))) #2 times length bc diploid
  snps.filt.cov.BD$D.REF.freq [i] <- (snps.filt.cov.BD$D.REF.sum[i] / (2*length(D.pop))) #2 times length bc diploid
  
  #alts
  
  snps.filt.cov$A.B.suc.ALT.freq[i] <- 1 - snps.filt.cov$A.B.suc.REF.freq[i]
  snps.filt.cov$C.D.res.ALT.freq[i] <- 1 - snps.filt.cov$C.D.res.REF.freq[i]
 
  snps.filt.cov.AC$A.ALT.freq[i] <- 1 - snps.filt.cov.AC$A.REF.freq[i] 
  snps.filt.cov.AC$C.ALT.freq [i] <- 1 - snps.filt.cov.AC$C.REF.freq[i]
  
  snps.filt.cov.AD$A.ALT.freq[i] <- 1 - snps.filt.cov.AD$A.REF.freq[i] 
  snps.filt.cov.AD$D.ALT.freq [i] <- 1 - snps.filt.cov.AD$D.REF.freq[i] 
  
  snps.filt.cov.BC$B.ALT.freq[i] <- 1 - snps.filt.cov.BC$B.REF.freq[i] 
  snps.filt.cov.BC$C.ALT.freq [i] <- 1 - snps.filt.cov.BC$C.REF.freq[i] 
  
  snps.filt.cov.BD$B.ALT.freq[i] <- 1 - snps.filt.cov.BD$B.REF.freq[i] 
  snps.filt.cov.BD$D.ALT.freq [i] <- 1 - snps.filt.cov.BD$D.REF.freq[i] 
  
  
  }


snps.filt.cov.AC$A.vs.C.REF.freq <- 0
snps.filt.cov.BC$B.vs.C.REF.freq <- 0
snps.filt.cov.AD$A.vs.D.REF.freq <- 0
snps.filt.cov.BD$B.vs.D.REF.freq <- 0

snps.filt.cov$Suc.AB.vs.Res.CD.REF.freq <- 0

#calculate abs difference in allele freqs between populations
for (i in 1:length(snps.filt.cov$GENE)) {
  snps.filt.cov.AC$A.vs.C.REF.freq[i] <- abs(as.numeric(snps.filt.cov.AC$A.REF.freq[i]) - as.numeric(snps.filt.cov.AC$C.REF.freq[i])) 
  snps.filt.cov.BC$B.vs.C.REF.freq[i] <- abs(as.numeric(snps.filt.cov.BC$B.REF.freq[i]) - as.numeric(snps.filt.cov.BC$C.REF.freq[i])) 
  snps.filt.cov.AD$A.vs.D.REF.freq[i] <- abs(as.numeric(snps.filt.cov.AD$A.REF.freq[i]) - as.numeric(snps.filt.cov.AD$D.REF.freq[i])) 
  snps.filt.cov.BD$B.vs.D.REF.freq[i] <- abs(as.numeric(snps.filt.cov.BD$B.REF.freq[i]) - as.numeric(snps.filt.cov.BD$D.REF.freq[i])) 
  
  snps.filt.cov$Suc.AB.vs.Res.CD.REF.freq[i] <- abs(as.numeric(snps.filt.cov$A.B.suc.REF.freq[i]) - as.numeric(snps.filt.cov$C.D.res.REF.freq[i])) 
  
  }


#filter out alleles where there is no difference
snps.filt.cov.var.AvC <- filter(snps.filt.cov.AC, A.vs.C.REF.freq != 0)
snps.filt.cov.var.BvC <- filter(snps.filt.cov.BC, B.vs.C.REF.freq != 0)
snps.filt.cov.var.AvD <- filter(snps.filt.cov.AD, A.vs.D.REF.freq != 0)
snps.filt.cov.var.BvD <- filter(snps.filt.cov.BD, B.vs.D.REF.freq != 0)

snps.filt.cov.var.SucvRes <- filter(snps.filt.cov, Suc.AB.vs.Res.CD.REF.freq != 0)


#order by largest to smallest difference
snps.filt.cov.var.AvC.sort <- snps.filt.cov.var.AvC[order(-snps.filt.cov.var.AvC$A.vs.C.REF.freq),]
snps.filt.cov.var.BvC.sort <- snps.filt.cov.var.BvC[order(-snps.filt.cov.var.BvC$B.vs.C.REF.freq),]
snps.filt.cov.var.AvD.sort <- snps.filt.cov.var.AvD[order(-snps.filt.cov.var.AvD$A.vs.D.REF.freq),]
snps.filt.cov.var.BvD.sort <- snps.filt.cov.var.BvD[order(-snps.filt.cov.var.BvD$B.vs.D.REF.freq),]

snps.filt.cov.var.SucvRes.sort <- snps.filt.cov.var.SucvRes[order(-snps.filt.cov.var.SucvRes$Suc.AB.vs.Res.CD.REF.freq),]



#Visualize difference in freqs
hist(snps.filt.cov.var.AvC.sort$A.vs.C.REF.freq, xlab = "Absolute difference between SNP frequencies for reference allele", main = "Population A (Tulare, suc) vs. Population C (General Beale, res)")

#visualize by base
ggplot(snps.filt.cov.var.AvC.sort, aes(x=A.vs.C.REF.freq)) +geom_histogram() + facet_grid(~REF) +xlab("Absolute difference between SNP frequencies for reference allele") + ggtitle("Population A (Tulare, suc) vs. Population C (General Beale, res)")
#ggsave(filename = 'plots/A.Tulare.Suc.vs.C.GeneralBeale.Res.Ref.Freq.pdf', plot = last_plot(), device = 'pdf', width = 8, height = 4, dpi = 300)


#Visualize individual freqs & by base
ggplot(snps.filt.cov.var.AvC.sort, aes(x=A.REF.freq)) +geom_histogram() + facet_grid(~REF) +xlab("SNP frequences for reference allele") + ggtitle("Population A (Tulare, suc)")
#ggsave(filename = 'plots/A.Tulare.Suc.Ref.Freq.pdf', plot = last_plot(), device = 'pdf', width = 8, height = 4, dpi = 300)


ggplot(snps.filt.cov.var.AvC.sort, aes(x=C.REF.freq)) +geom_histogram() + facet_grid(~REF) +xlab("SNP frequences for reference allele") + ggtitle("Population C (General Beale, res)")
#ggsave(filename = 'plots/C.GeneralBeale.Res.Ref.Freq.pdf', plot = last_plot(), device = 'pdf', width = 8, height = 4, dpi = 300)


#save file as csv
#write.csv(snps.filt.cov.var.AvC.sort, file = "results/GWSS_RNASeq1.snpEff.matrix_high.AandC_R_AlleleFreqs.csv")
#write.csv(snps.filt.cov.var.BvC.sort, file = "results/GWSS_RNASeq1.snpEff.matrix_high.BandC_R_AlleleFreqs.csv")
#write.csv(snps.filt.cov.var.AvD.sort, file = "results/GWSS_RNASeq1.snpEff.matrix_high.AandD_R_AlleleFreqs.csv")
#write.csv(snps.filt.cov.var.BvD.sort, file = "results/GWSS_RNASeq1.snpEff.matrix_high.BandD_R_AlleleFreqs.csv")
#write.csv(snps.filt.cov.var.SucvRes.sort, file = "results/GWSS_RNASeq1.snpEff.matrix_high.SucvRes_R_AlleleFreqs.csv")


#make short file
#remove some columns generated in R
snps.filt.cov.var.AvC.sort.short <- subset(snps.filt.cov.var.AvC.sort, select = -c(A.REF.sum, C.REF.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop))
snps.filt.cov.var.BvC.sort.short <- subset(snps.filt.cov.var.BvC.sort, select = -c(B.REF.sum, C.REF.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop))
snps.filt.cov.var.BvD.sort.short <- subset(snps.filt.cov.var.BvD.sort, select = -c(B.REF.sum, D.REF.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop))
snps.filt.cov.var.AvD.sort.short <- subset(snps.filt.cov.var.AvD.sort, select = -c(A.REF.sum, D.REF.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop))
snps.filt.cov.var.SucvRes.sort.short <- subset(snps.filt.cov.var.SucvRes.sort, select = -c(A.REF.sum, B.REF.sum, D.REF.sum, C.REF.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop))


#filter out alleles that don't differ > .5
snps.filt.cov.var.AvC.sort.short.filt <- filter(snps.filt.cov.var.AvC.sort.short, A.vs.C.REF.freq > 0.5)
snps.filt.cov.var.BvC.sort.short.filt <- filter(snps.filt.cov.var.BvC.sort.short, B.vs.C.REF.freq > 0.5)
snps.filt.cov.var.BvD.sort.short.filt <- filter(snps.filt.cov.var.BvD.sort.short, B.vs.D.REF.freq > 0.5)
snps.filt.cov.var.AvD.sort.short.filt <- filter(snps.filt.cov.var.AvD.sort.short, A.vs.D.REF.freq > 0.5)
snps.filt.cov.var.SucvRes.sort.short.filt <- filter(snps.filt.cov.var.SucvRes.sort.short, Suc.AB.vs.Res.CD.REF.freq > 0.5)


#save file as csv
# write.csv(snps.filt.cov.var.AvC.sort.short.filt, file = "results/GWSS_RNASeq1.snpEff.matrix_high.AandC_R_AlleleFreqs.short.csv")
# write.csv(snps.filt.cov.var.BvC.sort.short.filt, file = "results/GWSS_RNASeq1.snpEff.matrix_high.BandC_R_AlleleFreqs.short.csv")
# write.csv(snps.filt.cov.var.BvD.sort.short.filt, file = "results/GWSS_RNASeq1.snpEff.matrix_high.BandD_R_AlleleFreqs.short.csv")
# write.csv(snps.filt.cov.var.AvD.sort.short.filt, file = "results/GWSS_RNASeq1.snpEff.matrix_high.AandD_R_AlleleFreqs.short.csv")
# write.csv(snps.filt.cov.var.SucvRes.sort.short.filt, file = "results/GWSS_RNASeq1.snpEff.matrix_high.SucandRes_R_AlleleFreqs.short.csv")



#cross reference the different lists for overlap
a <- full_join(snps.filt.cov.var.SucvRes.sort.short, snps.filt.cov.var.BvC.sort.short)
b <- full_join(a, snps.filt.cov.var.BvD.sort.short)
c <- full_join(b, snps.filt.cov.var.AvD.sort.short)
d <- full_join(c, snps.filt.cov.var.AvC.sort.short)

#rename
d$FreqDiff_Susceptible_AB_vs_Resistant_CD <- d$Suc.AB.vs.Res.CD.REF.freq
d$FreqDiff_Tulare_Susceptible_A_vs_GeneralBeale_Resistant_C <- d$A.vs.C.REF.freq
d$FreqDiff_Tulare_Susceptible_A_vs_Tulare_Resistant_D <- d$A.vs.D.REF.freq
d$FreqDiff_Temecula_Susceptible_B_vs_GeneralBeale_Resistant_C <- d$B.vs.C.REF.freq
d$FreqDiff_Temecula_Susceptible_B_vs_Tulare_Resistant_D <- d$B.vs.D.REF.freq


#subset to only get differences between reference freq for group comparions
d.v2 <- subset(d, select = -c(A.REF.freq, B.REF.freq, C.REF.freq, D.REF.freq, A.ALT.freq, B.ALT.freq, C.ALT.freq,D.ALT.freq,A.B.suc.ALT.freq, A.B.suc.REF.freq, C.D.res.ALT.freq, C.D.res.REF.freq, Suc.AB.vs.Res.CD.REF.freq, A.vs.C.REF.freq, A.vs.D.REF.freq, B.vs.D.REF.freq, B.vs.C.REF.freq))

#reorder
d.v2 <- d.v2[order(-d.v2$FreqDiff_Susceptible_AB_vs_Resistant_CD),]


#write.csv(d.v2, file = "results/GWSS_RNASeq1.snpEff.matrix_high.CombinedResults_R_AlleleFreqsDifferences.csv")


#filter
d.v2.filt <- filter(d.v2, FreqDiff_Susceptible_AB_vs_Resistant_CD > 0.5 | FreqDiff_Tulare_Susceptible_A_vs_GeneralBeale_Resistant_C > 0.625 | FreqDiff_Tulare_Susceptible_A_vs_Tulare_Resistant_D > 0.625 | FreqDiff_Temecula_Susceptible_B_vs_GeneralBeale_Resistant_C > 0.625 | FreqDiff_Temecula_Susceptible_B_vs_Tulare_Resistant_D > 0.625)


#write.csv(d.v2.filt, file = "results/GWSS_RNASeq1.snpEff.matrix_high.CombinedResults_R_AlleleFreqsDifferences.short.csv")
