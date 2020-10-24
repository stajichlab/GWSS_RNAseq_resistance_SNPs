#GWSS Indel analysis between resistant & susceptible populations
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
indel.filt <- filter(snps, indel == TRUE)


#filter to only get indels sequenced across all individuals 
#I expect this will probably filter out many indels 
#but with only 4 reps per population, trying to be conservative

#set up column with all data assumed to be present
indel.filt$missing <- FALSE
indel.filt$missing.Apop <- FALSE
indel.filt$missing.Bpop <- FALSE
indel.filt$missing.Cpop <- FALSE
indel.filt$missing.Dpop <- FALSE

#for each individual, check if indel called, if missing - reassign value to TRUE
#else leave value as is, with default of FALSE
for (pop in pop_col) {
  
  indel.filt %<>%
    mutate(missing = ifelse(indel.filt[pop] == ".", TRUE, as.character(missing)))
  
  if (pop %in% A.pop) {
    indel.filt %<>%
      mutate(missing.Apop = ifelse(indel.filt[pop] == ".", TRUE, as.character(missing.Apop)))
  } else if (pop %in% B.pop) {
    indel.filt %<>%
      mutate(missing.Bpop = ifelse(indel.filt[pop] == ".", TRUE, as.character(missing.Bpop)))
    
  } else if (pop %in% C.pop) {
    indel.filt %<>%
      mutate(missing.Cpop = ifelse(indel.filt[pop] == ".", TRUE, as.character(missing.Cpop)))
    
  } else {
    indel.filt %<>%
      mutate(missing.Dpop = ifelse(indel.filt[pop] == ".", TRUE, as.character(missing.Dpop)))
    
  }
}


#assessing number of *s
indel.filt$notcalled <- FALSE

for (pop in pop_col) {
  for (i in 1:length(indel.filt$GENE)) {
    if (str_count(indel.filt[,pop][i], pattern = "\\*") > 0) {
      indel.filt$notcalled[i] <- TRUE
    } } }



#assessing if more than one alternate allele
indel.filt$alt.multalts <- FALSE

for (i in 1:length(indel.filt$GENE)) {
  if (str_count(indel.filt$ALT[i], pattern = ",") > 0) {
      indel.filt$alt.multalts[i] <- TRUE
     } }


#filter out all rows with *s and multiple alternate alleles 
indel.filt <- filter(indel.filt, notcalled == FALSE & alt.multalts == FALSE)


#we need to always sum with the longer allele 
#this is to help avoid overcounting
#though we may still run into issues with this with double counting
#lets see!
indel.filt$usealtsum <- FALSE

for (i in 1:length(indel.filt$GENE)) {
  if (nchar(indel.filt$REF[i]) < nchar(indel.filt$ALT[i])) {
    indel.filt$usealtsum[i] <- TRUE
  }
}


#filter to only have table where indels were covered for all individuals across all pops
indel.filt.cov <- filter(indel.filt, missing == FALSE)

#filter to only have table where SNPS were covered for all individuals in comparison being made
indel.filt.cov.AC <- filter(indel.filt, missing.Apop == FALSE &  missing.Cpop == FALSE)
indel.filt.cov.AD <- filter(indel.filt, missing.Apop == FALSE &  missing.Dpop == FALSE)
indel.filt.cov.BC <- filter(indel.filt, missing.Bpop == FALSE &  missing.Cpop == FALSE)
indel.filt.cov.BD <- filter(indel.filt, missing.Bpop == FALSE &  missing.Dpop == FALSE)


# /:  genotype unphased◦ |:  genotype phased
# I think we can ignore these for now as we only want allele freqs

#calculate allele frequencies for reference allele 
#we need to deal with *s

#set up columns to save sums for whole dataset
indel.filt.cov$A.REF.sum <- 0
indel.filt.cov$B.REF.sum <- 0
indel.filt.cov$C.REF.sum <- 0
indel.filt.cov$D.REF.sum <- 0

#set up columns to save sums for each pop comparison
indel.filt.cov.AC$A.REF.sum <- 0
indel.filt.cov.AC$C.REF.sum <- 0

indel.filt.cov.AD$A.REF.sum <- 0
indel.filt.cov.AD$D.REF.sum <- 0

indel.filt.cov.BC$B.REF.sum <- 0
indel.filt.cov.BC$C.REF.sum <- 0

indel.filt.cov.BD$B.REF.sum <- 0
indel.filt.cov.BD$D.REF.sum <- 0

#set up alt sums
indel.filt.cov$A.ALT.sum <- 0
indel.filt.cov$B.ALT.sum <- 0
indel.filt.cov$C.ALT.sum <- 0
indel.filt.cov$D.ALT.sum <- 0

indel.filt.cov.AC$A.ALT.sum <- 0
indel.filt.cov.AC$C.ALT.sum <- 0

indel.filt.cov.AD$A.ALT.sum <- 0
indel.filt.cov.AD$D.ALT.sum <- 0

indel.filt.cov.BC$B.ALT.sum <- 0
indel.filt.cov.BC$C.ALT.sum <- 0

indel.filt.cov.BD$B.ALT.sum <- 0
indel.filt.cov.BD$D.ALT.sum <- 0


  
#sum A pop ref alleles
for (pop in A.pop) {
  for (i in 1:length(indel.filt.cov$GENE)) {
    if (indel.filt.cov$usealtsum[i] == FALSE) {  
    indel.filt.cov$A.REF.sum[i] <- indel.filt.cov$A.REF.sum[i] + str_count(indel.filt.cov[,pop][i], pattern = indel.filt.cov$REF[i])
    }
    else {
    indel.filt.cov$A.ALT.sum[i] <- indel.filt.cov$A.ALT.sum[i] + str_count(indel.filt.cov[,pop][i], pattern = indel.filt.cov$ALT[i])
    } }
  for (i in 1:length(indel.filt.cov.AC$GENE)) {
    if (indel.filt.cov.AC$usealtsum[i] == FALSE) {  
     indel.filt.cov.AC$A.REF.sum[i] <- indel.filt.cov.AC$A.REF.sum[i] + str_count(indel.filt.cov.AC[,pop][i], pattern = indel.filt.cov.AC$REF[i])
    }
    else {
     indel.filt.cov.AC$A.ALT.sum[i] <- indel.filt.cov.AC$A.ALT.sum[i] + str_count(indel.filt.cov.AC[,pop][i], pattern = indel.filt.cov.AC$ALT[i])
    }
  }
  for (i in 1:length(indel.filt.cov.AD$GENE)) {
    if (indel.filt.cov.AD$usealtsum[i] == FALSE) {  
    indel.filt.cov.AD$A.REF.sum[i] <- indel.filt.cov.AD$A.REF.sum[i] + str_count(indel.filt.cov.AD[,pop][i], pattern = indel.filt.cov.AD$REF[i])
    } else {
    indel.filt.cov.AD$A.ALT.sum[i] <- indel.filt.cov.AD$A.ALT.sum[i] + str_count(indel.filt.cov.AD[,pop][i], pattern = indel.filt.cov.AD$ALT[i])
    }
      }
}

#sum B pop ref alleles
for (pop in B.pop) {
  for (i in 1:length(indel.filt.cov$GENE)) {
    if (indel.filt.cov$usealtsum[i] == FALSE) {  
    indel.filt.cov$B.REF.sum[i] <- indel.filt.cov$B.REF.sum[i] + str_count(indel.filt.cov[,pop][i], pattern = indel.filt.cov$REF[i])
    } else {
      indel.filt.cov$B.ALT.sum[i] <- indel.filt.cov$B.ALT.sum[i] + str_count(indel.filt.cov[,pop][i], pattern = indel.filt.cov$ALT[i])
      
    }
    }
  for (i in 1:length(indel.filt.cov.BC$GENE)) {
    if (indel.filt.cov.BC$usealtsum[i] == FALSE) {  
    indel.filt.cov.BC$B.REF.sum[i] <- indel.filt.cov.BC$B.REF.sum[i] + str_count(indel.filt.cov.BC[,pop][i], pattern = indel.filt.cov.BC$REF[i])
    } else {
      indel.filt.cov.BC$B.ALT.sum[i] <- indel.filt.cov.BC$B.ALT.sum[i] + str_count(indel.filt.cov.BC[,pop][i], pattern = indel.filt.cov.BC$ALT[i])
      
    }
     }
  for (i in 1:length(indel.filt.cov.BD$GENE)) {
    if (indel.filt.cov.BD$usealtsum[i] == FALSE) {  
      indel.filt.cov.BD$B.REF.sum[i] <- indel.filt.cov.BD$B.REF.sum[i] + str_count(indel.filt.cov.BD[,pop][i], pattern = indel.filt.cov.BD$REF[i])
    } else {
      indel.filt.cov.BD$B.ALT.sum[i] <- indel.filt.cov.BD$B.ALT.sum[i] + str_count(indel.filt.cov.BD[,pop][i], pattern = indel.filt.cov.BD$ALT[i])
      
  } }
}

#sum C pop ref alleles
for (pop in C.pop) {
  for (i in 1:length(indel.filt.cov$GENE)) {
    if (indel.filt.cov$usealtsum[i] == FALSE) {  
    indel.filt.cov$C.REF.sum[i] <- indel.filt.cov$C.REF.sum[i] + str_count(indel.filt.cov[,pop][i], pattern = indel.filt.cov$REF[i])
    } else {
    indel.filt.cov$C.ALT.sum[i] <- indel.filt.cov$C.ALT.sum[i] + str_count(indel.filt.cov[,pop][i], pattern = indel.filt.cov$ALT[i])
      
  }}
  for (i in 1:length(indel.filt.cov.BC$GENE)) {
    if (indel.filt.cov.BC$usealtsum[i] == FALSE) {  
    indel.filt.cov.BC$C.REF.sum[i] <- indel.filt.cov.BC$C.REF.sum[i] + str_count(indel.filt.cov.BC[,pop][i], pattern = indel.filt.cov.BC$REF[i])
    } else {
      indel.filt.cov.BC$C.ALT.sum[i] <- indel.filt.cov.BC$C.ALT.sum[i] + str_count(indel.filt.cov.BC[,pop][i], pattern = indel.filt.cov.BC$ALT[i])
    
  }}
  for (i in 1:length(indel.filt.cov.AC$GENE)) {
    if (indel.filt.cov.AC$usealtsum[i] == FALSE) {  
    indel.filt.cov.AC$C.REF.sum[i] <- indel.filt.cov.AC$C.REF.sum[i] + str_count(indel.filt.cov.AC[,pop][i], pattern = indel.filt.cov.AC$REF[i])
    } else {
    indel.filt.cov.AC$C.ALT.sum[i] <- indel.filt.cov.AC$C.ALT.sum[i] + str_count(indel.filt.cov.AC[,pop][i], pattern = indel.filt.cov.AC$ALT[i])
      
  }}
}


#sum D pop ref alleles
for (pop in D.pop) {
  for (i in 1:length(indel.filt.cov$GENE)) {
    if (indel.filt.cov$usealtsum[i] == FALSE) {  
    indel.filt.cov$D.REF.sum[i] <- indel.filt.cov$D.REF.sum[i] + str_count(indel.filt.cov[,pop][i], pattern = indel.filt.cov$REF[i])
    } else {
    indel.filt.cov$D.ALT.sum[i] <- indel.filt.cov$D.ALT.sum[i] + str_count(indel.filt.cov[,pop][i], pattern = indel.filt.cov$ALT[i])
      
  }}
  for (i in 1:length(indel.filt.cov.BD$GENE)) {
    if (indel.filt.cov.BD$usealtsum[i] == FALSE) {  
    indel.filt.cov.BD$D.REF.sum[i] <- indel.filt.cov.BD$D.REF.sum[i] + str_count(indel.filt.cov.BD[,pop][i], pattern = indel.filt.cov.BD$REF[i])
    } else {
      indel.filt.cov.BD$D.ALT.sum[i] <- indel.filt.cov.BD$D.ALT.sum[i] + str_count(indel.filt.cov.BD[,pop][i], pattern = indel.filt.cov.BD$ALT[i])
      
  }}
  for (i in 1:length(indel.filt.cov.AD$GENE)) {
    if (indel.filt.cov.AD$usealtsum[i] == FALSE) {  
    indel.filt.cov.AD$D.REF.sum[i] <- indel.filt.cov.AD$D.REF.sum[i] + str_count(indel.filt.cov.AD[,pop][i], pattern = indel.filt.cov.AD$REF[i])
    } else {
      indel.filt.cov.AD$D.ALT.sum[i] <- indel.filt.cov.AD$D.ALT.sum[i] + str_count(indel.filt.cov.AD[,pop][i], pattern = indel.filt.cov.AD$ALT[i])
      
  }}
}

#A - Tulare (suc)
#B - Temecula (suc)
#C - GeneralBeale (res)
#D - Tulare (res)


#set up to save lumped by resistance status
indel.filt.cov$A.B.suc.REF.freq <- 0
indel.filt.cov$C.D.res.REF.freq <-0

indel.filt.cov$A.B.suc.ALT.freq <- 0
indel.filt.cov$C.D.res.ALT.freq <-0


#set up to save freq for each pop 
indel.filt.cov.AC$A.REF.freq <- 0
indel.filt.cov.AC$C.REF.freq <- 0

indel.filt.cov.AD$A.REF.freq <- 0
indel.filt.cov.AD$D.REF.freq <- 0

indel.filt.cov.BC$B.REF.freq <- 0
indel.filt.cov.BC$C.REF.freq <- 0

indel.filt.cov.BD$B.REF.freq <- 0
indel.filt.cov.BD$D.REF.freq <- 0


#alt freqs
indel.filt.cov.AC$A.ALT.freq <- 0
indel.filt.cov.AC$C.ALT.freq <- 0

indel.filt.cov.AD$A.ALT.freq <- 0
indel.filt.cov.AD$D.ALT.freq <- 0

indel.filt.cov.BC$B.ALT.freq <- 0
indel.filt.cov.BC$C.ALT.freq <- 0

indel.filt.cov.BD$B.ALT.freq <- 0
indel.filt.cov.BD$D.ALT.freq <- 0





#turn sums to proportion & calculate alt allele freq based on REF
for (i in 1:length(indel.filt.cov$GENE)) {
  
  if (indel.filt.cov$usealtsum[i] == FALSE) {
    indel.filt.cov$A.B.suc.REF.freq[i] <- ((indel.filt.cov$A.REF.sum[i]+indel.filt.cov$B.REF.sum[i]) / (2*(length(A.pop)+length(B.pop)))) 
    indel.filt.cov$C.D.res.REF.freq[i] <- ((indel.filt.cov$C.REF.sum[i]+indel.filt.cov$D.REF.sum[i]) / (2*(length(C.pop)+length(D.pop)))) 
    
    indel.filt.cov$A.B.suc.ALT.freq[i] <- 1 - indel.filt.cov$A.B.suc.REF.freq[i]
    indel.filt.cov$C.D.res.ALT.freq[i] <- 1 - indel.filt.cov$C.D.res.REF.freq[i]
  } else {
    indel.filt.cov$A.B.suc.ALT.freq[i] <- ((indel.filt.cov$A.ALT.sum[i]+indel.filt.cov$B.ALT.sum[i]) / (2*(length(A.pop)+length(B.pop)))) 
    indel.filt.cov$C.D.res.ALT.freq[i] <- ((indel.filt.cov$C.ALT.sum[i]+indel.filt.cov$D.ALT.sum[i]) / (2*(length(C.pop)+length(D.pop)))) 
    
    indel.filt.cov$A.B.suc.REF.freq[i] <- 1 - indel.filt.cov$A.B.suc.ALT.freq[i]
    indel.filt.cov$C.D.res.REF.freq[i] <- 1 - indel.filt.cov$C.D.res.ALT.freq[i]
    
  }
  
  if (indel.filt.cov.AC$usealtsum[i] == FALSE) {
    indel.filt.cov.AC$A.REF.freq[i] <- (indel.filt.cov.AC$A.REF.sum[i] / (2*length(A.pop))) #2 times length bc diploid
    indel.filt.cov.AC$C.REF.freq [i] <- (indel.filt.cov.AC$C.REF.sum[i] / (2*length(C.pop))) #2 times length bc diploid
    
    indel.filt.cov.AC$A.ALT.freq[i] <- 1 - indel.filt.cov.AC$A.REF.freq[i] 
    indel.filt.cov.AC$C.ALT.freq [i] <- 1 - indel.filt.cov.AC$C.REF.freq[i]
    
     } else {
    indel.filt.cov.AC$A.ALT.freq[i] <- (indel.filt.cov.AC$A.ALT.sum[i] / (2*length(A.pop))) #2 times length bc diploid
    indel.filt.cov.AC$C.ALT.freq [i] <- (indel.filt.cov.AC$C.ALT.sum[i] / (2*length(C.pop))) #2 times length bc diploid
    
    indel.filt.cov.AC$A.REF.freq[i] <- 1 - indel.filt.cov.AC$A.ALT.freq[i] 
    indel.filt.cov.AC$C.REF.freq [i] <- 1 - indel.filt.cov.AC$C.ALT.freq[i]
    
  }
  
  if (indel.filt.cov.AD$usealtsum[i] == FALSE) {
    indel.filt.cov.AD$A.REF.freq[i] <- (indel.filt.cov.AD$A.REF.sum[i] / (2*length(A.pop))) #2 times length bc diploid
    indel.filt.cov.AD$D.REF.freq [i] <- (indel.filt.cov.AD$D.REF.sum[i] / (2*length(D.pop))) #2 times length bc diploid
    
    indel.filt.cov.AD$A.ALT.freq[i] <- 1 - indel.filt.cov.AD$A.REF.freq[i] 
    indel.filt.cov.AD$D.ALT.freq [i] <- 1 - indel.filt.cov.AD$D.REF.freq[i] 
    
  } else {
    indel.filt.cov.AD$A.ALT.freq[i] <- (indel.filt.cov.AD$A.ALT.sum[i] / (2*length(A.pop))) #2 times length bc diploid
    indel.filt.cov.AD$D.ALT.freq [i] <- (indel.filt.cov.AD$D.ALT.sum[i] / (2*length(D.pop))) #2 times length bc diploid
    
    indel.filt.cov.AD$A.REF.freq[i] <- 1 - indel.filt.cov.AD$A.ALT.freq[i] 
    indel.filt.cov.AD$D.REF.freq [i] <- 1 - indel.filt.cov.AD$D.ALT.freq[i] 
    
  }
 
  if (indel.filt.cov.BC$usealtsum[i] == FALSE) {
    indel.filt.cov.BC$B.REF.freq[i] <- (indel.filt.cov.BC$B.REF.sum[i] / (2*length(B.pop))) #2 times length bc diploid
    indel.filt.cov.BC$C.REF.freq [i] <- (indel.filt.cov.BC$C.REF.sum[i] / (2*length(C.pop))) #2 times length bc diploid
    
    indel.filt.cov.BC$B.ALT.freq[i] <- 1 - indel.filt.cov.BC$B.REF.freq[i] 
    indel.filt.cov.BC$C.ALT.freq [i] <- 1 - indel.filt.cov.BC$C.REF.freq[i] 
    
  } else {
    indel.filt.cov.BC$B.ALT.freq[i] <- (indel.filt.cov.BC$B.ALT.sum[i] / (2*length(B.pop))) #2 times length bc diploid
    indel.filt.cov.BC$C.ALT.freq [i] <- (indel.filt.cov.BC$C.ALT.sum[i] / (2*length(C.pop))) #2 times length bc diploid
    
    indel.filt.cov.BC$B.REF.freq[i] <- 1 - indel.filt.cov.BC$B.ALT.freq[i] 
    indel.filt.cov.BC$C.REF.freq [i] <- 1 - indel.filt.cov.BC$C.ALT.freq[i] 
    
  }

  if (indel.filt.cov.BD$usealtsum[i] == FALSE) {
    indel.filt.cov.BD$B.REF.freq[i] <- (indel.filt.cov.BD$B.REF.sum[i] / (2*length(B.pop))) #2 times length bc diploid
    indel.filt.cov.BD$D.REF.freq [i] <- (indel.filt.cov.BD$D.REF.sum[i] / (2*length(D.pop))) #2 times length bc diploid
    
    indel.filt.cov.BD$B.ALT.freq[i] <- 1 - indel.filt.cov.BD$B.REF.freq[i] 
    indel.filt.cov.BD$D.ALT.freq [i] <- 1 - indel.filt.cov.BD$D.REF.freq[i] 
    
  } else {
    indel.filt.cov.BD$B.ALT.freq[i] <- (indel.filt.cov.BD$B.ALT.sum[i] / (2*length(B.pop))) #2 times length bc diploid
    indel.filt.cov.BD$D.ALT.freq [i] <- (indel.filt.cov.BD$D.ALT.sum[i] / (2*length(D.pop))) #2 times length bc diploid
    
    indel.filt.cov.BD$B.REF.freq[i] <- 1 - indel.filt.cov.BD$B.ALT.freq[i] 
    indel.filt.cov.BD$D.REF.freq [i] <- 1 - indel.filt.cov.BD$D.ALT.freq[i] 
    
  }
  
}


indel.filt.cov.AC$A.vs.C.REF.freq <- 0
indel.filt.cov.BC$B.vs.C.REF.freq <- 0
indel.filt.cov.AD$A.vs.D.REF.freq <- 0
indel.filt.cov.BD$B.vs.D.REF.freq <- 0

indel.filt.cov$Suc.AB.vs.Res.CD.REF.freq <- 0

#calculate abs difference in allele freqs between populations
for (i in 1:length(indel.filt.cov$GENE)) {
  indel.filt.cov.AC$A.vs.C.REF.freq[i] <- abs(as.numeric(indel.filt.cov.AC$A.REF.freq[i]) - as.numeric(indel.filt.cov.AC$C.REF.freq[i])) 
  indel.filt.cov.BC$B.vs.C.REF.freq[i] <- abs(as.numeric(indel.filt.cov.BC$B.REF.freq[i]) - as.numeric(indel.filt.cov.BC$C.REF.freq[i])) 
  indel.filt.cov.AD$A.vs.D.REF.freq[i] <- abs(as.numeric(indel.filt.cov.AD$A.REF.freq[i]) - as.numeric(indel.filt.cov.AD$D.REF.freq[i])) 
  indel.filt.cov.BD$B.vs.D.REF.freq[i] <- abs(as.numeric(indel.filt.cov.BD$B.REF.freq[i]) - as.numeric(indel.filt.cov.BD$D.REF.freq[i])) 
  
  indel.filt.cov$Suc.AB.vs.Res.CD.REF.freq[i] <- abs(as.numeric(indel.filt.cov$A.B.suc.REF.freq[i]) - as.numeric(indel.filt.cov$C.D.res.REF.freq[i])) 
  
}


#filter out alleles where there is no difference
indel.filt.cov.var.AvC <- filter(indel.filt.cov.AC, A.vs.C.REF.freq != 0)
indel.filt.cov.var.BvC <- filter(indel.filt.cov.BC, B.vs.C.REF.freq != 0)
indel.filt.cov.var.AvD <- filter(indel.filt.cov.AD, A.vs.D.REF.freq != 0)
indel.filt.cov.var.BvD <- filter(indel.filt.cov.BD, B.vs.D.REF.freq != 0)

indel.filt.cov.var.SucvRes <- filter(indel.filt.cov, Suc.AB.vs.Res.CD.REF.freq != 0)


#order by largest to smallest difference
indel.filt.cov.var.AvC.sort <- indel.filt.cov.var.AvC[order(-indel.filt.cov.var.AvC$A.vs.C.REF.freq),]
indel.filt.cov.var.BvC.sort <- indel.filt.cov.var.BvC[order(-indel.filt.cov.var.BvC$B.vs.C.REF.freq),]
indel.filt.cov.var.AvD.sort <- indel.filt.cov.var.AvD[order(-indel.filt.cov.var.AvD$A.vs.D.REF.freq),]
indel.filt.cov.var.BvD.sort <- indel.filt.cov.var.BvD[order(-indel.filt.cov.var.BvD$B.vs.D.REF.freq),]

indel.filt.cov.var.SucvRes.sort <- indel.filt.cov.var.SucvRes[order(-indel.filt.cov.var.SucvRes$Suc.AB.vs.Res.CD.REF.freq),]



#Visualize difference in freqs
hist(indel.filt.cov.var.AvC.sort$A.vs.C.REF.freq, xlab = "Absolute difference between indel frequencies for reference allele", main = "Population A (Tulare, suc) vs. Population C (General Beale, res)")

#save file as csv
#write.csv(indel.filt.cov.var.AvC.sort, file = "intermed_results/GWSS_RNASeq1.snpEff.matrix_high.AandC_R_indelfreqs.csv")
#write.csv(indel.filt.cov.var.BvC.sort, file = "intermed_results/GWSS_RNASeq1.snpEff.matrix_high.BandC_R_indelfreqs.csv")
#write.csv(indel.filt.cov.var.AvD.sort, file = "intermed_results/GWSS_RNASeq1.snpEff.matrix_high.AandD_R_indelfreqs.csv")
#write.csv(indel.filt.cov.var.BvD.sort, file = "intermed_results/GWSS_RNASeq1.snpEff.matrix_high.BandD_R_indelfreqs.csv")
#write.csv(indel.filt.cov.var.SucvRes.sort, file = "intermed_results/GWSS_RNASeq1.snpEff.matrix_high.SucvRes_R_indelfreqs.csv")


#make short file
#remove some columns generated in R
indel.filt.cov.var.AvC.sort.short <- subset(indel.filt.cov.var.AvC.sort, select = -c(A.REF.sum, C.REF.sum,A.ALT.sum, C.ALT.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop, notcalled, alt.multalts, usealtsum))
indel.filt.cov.var.BvC.sort.short <- subset(indel.filt.cov.var.BvC.sort, select = -c(B.REF.sum, C.REF.sum,B.ALT.sum, C.ALT.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop, notcalled, alt.multalts, usealtsum))
indel.filt.cov.var.BvD.sort.short <- subset(indel.filt.cov.var.BvD.sort, select = -c(B.REF.sum, D.REF.sum, B.ALT.sum, D.ALT.sum,indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop, notcalled, alt.multalts, usealtsum))
indel.filt.cov.var.AvD.sort.short <- subset(indel.filt.cov.var.AvD.sort, select = -c(A.REF.sum, D.REF.sum, A.ALT.sum, D.ALT.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop, notcalled, alt.multalts, usealtsum))
indel.filt.cov.var.SucvRes.sort.short <- subset(indel.filt.cov.var.SucvRes.sort, select = -c(A.REF.sum, B.REF.sum, D.REF.sum, C.REF.sum,A.ALT.sum, B.ALT.sum, D.ALT.sum, C.ALT.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop, notcalled, alt.multalts, usealtsum))


#filter out alleles that don't differ > .5
indel.filt.cov.var.AvC.sort.short.filt <- filter(indel.filt.cov.var.AvC.sort.short, A.vs.C.REF.freq > 0.5)
indel.filt.cov.var.BvC.sort.short.filt <- filter(indel.filt.cov.var.BvC.sort.short, B.vs.C.REF.freq > 0.5)
indel.filt.cov.var.BvD.sort.short.filt <- filter(indel.filt.cov.var.BvD.sort.short, B.vs.D.REF.freq > 0.5)
indel.filt.cov.var.AvD.sort.short.filt <- filter(indel.filt.cov.var.AvD.sort.short, A.vs.D.REF.freq > 0.5)
indel.filt.cov.var.SucvRes.sort.short.filt <- filter(indel.filt.cov.var.SucvRes.sort.short, Suc.AB.vs.Res.CD.REF.freq > 0.5)


#save file as csv
#write.csv(indel.filt.cov.var.AvC.sort.short.filt, file = "intermed_results/GWSS_RNASeq1.snpEff.matrix_high.AandC_R_indelfreqs.short.csv")
#write.csv(indel.filt.cov.var.BvC.sort.short.filt, file = "intermed_results/GWSS_RNASeq1.snpEff.matrix_high.BandC_R_indelfreqs.short.csv")
#write.csv(indel.filt.cov.var.BvD.sort.short.filt, file = "intermed_results/GWSS_RNASeq1.snpEff.matrix_high.BandD_R_indelfreqs.short.csv")
#write.csv(indel.filt.cov.var.AvD.sort.short.filt, file = "intermed_results/GWSS_RNASeq1.snpEff.matrix_high.AandD_R_indelfreqs.short.csv") 
#write.csv(indel.filt.cov.var.SucvRes.sort.short.filt, file = "intermed_results/GWSS_RNASeq1.snpEff.matrix_high.SucandRes_R_indelfreqs.short.csv")


#remove some columns generated in R
indel.filt.cov.var.AvC.short <- subset(indel.filt.cov.var.AvC, select = -c(A.REF.sum, C.REF.sum, A.ALT.sum, C.ALT.sum,indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop,notcalled, alt.multalts, usealtsum))
indel.filt.cov.var.BvC.short <- subset(indel.filt.cov.var.BvC, select = -c(B.REF.sum, C.REF.sum, B.ALT.sum, C.ALT.sum,indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop,notcalled, alt.multalts, usealtsum))
indel.filt.cov.var.BvD.short <- subset(indel.filt.cov.var.BvD, select = -c(B.REF.sum, D.REF.sum, B.ALT.sum, D.ALT.sum,indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop,notcalled, alt.multalts, usealtsum))
indel.filt.cov.var.AvD.short <- subset(indel.filt.cov.var.AvD, select = -c(A.REF.sum, D.REF.sum,A.ALT.sum,D.ALT.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop,notcalled, alt.multalts, usealtsum))
indel.filt.cov.var.SucvRes.short <- subset(indel.filt.cov.var.SucvRes, select = -c(A.REF.sum, B.REF.sum, D.REF.sum, C.REF.sum,A.ALT.sum, B.ALT.sum, D.ALT.sum, C.ALT.sum, indel, missing, missing.Apop, missing.Bpop, missing.Cpop, missing.Dpop,notcalled, alt.multalts, usealtsum))




#cross reference the different lists for overlap
a <- full_join(indel.filt.cov.var.SucvRes.short, indel.filt.cov.var.BvC.short)
b <- full_join(a, indel.filt.cov.var.BvD.short, by = c("CHROM", "POS", "FLANKING", "TYPE", "IMPACT", "GENE", "CHANGEDNA", "CHANGEPEP", "REF", "ALT", "A.1", "A.2", "A.3", "A.4", "B.1", "B.2", "B.3", "B.4", "C.1", "C.2", "C.3", "C.4", "D.1", "D.2", "D.3", "D.4", "ANN"))
c <- full_join(b, indel.filt.cov.var.AvD.short, by = c("CHROM", "POS", "FLANKING", "TYPE", "IMPACT", "GENE", "CHANGEDNA", "CHANGEPEP", "REF", "ALT", "A.1", "A.2", "A.3", "A.4", "B.1", "B.2", "B.3", "B.4", "C.1", "C.2", "C.3", "C.4", "D.1", "D.2", "D.3", "D.4", "ANN"))
d <- full_join(c, indel.filt.cov.var.AvC.short, by = c("CHROM", "POS", "FLANKING", "TYPE", "IMPACT", "GENE", "CHANGEDNA", "CHANGEPEP", "REF", "ALT", "A.1", "A.2", "A.3", "A.4", "B.1", "B.2", "B.3", "B.4", "C.1", "C.2", "C.3", "C.4", "D.1", "D.2", "D.3", "D.4", "ANN"))

#rename
d$FreqDiff_Susceptible_AB_vs_Resistant_CD <- d$Suc.AB.vs.Res.CD.REF.freq
d$FreqDiff_Tulare_Susceptible_A_vs_GeneralBeale_Resistant_C <- d$A.vs.C.REF.freq
d$FreqDiff_Tulare_Susceptible_A_vs_Tulare_Resistant_D <- d$A.vs.D.REF.freq
d$FreqDiff_Temecula_Susceptible_B_vs_GeneralBeale_Resistant_C <- d$B.vs.C.REF.freq
d$FreqDiff_Temecula_Susceptible_B_vs_Tulare_Resistant_D <- d$B.vs.D.REF.freq


#subset to only get differences between reference freq for group comparions
d.v2 <- subset(d, select = -c(A.REF.freq.x, A.REF.freq.y, B.REF.freq.x, B.REF.freq.y, C.REF.freq.x, C.REF.freq.y, D.REF.freq.x, D.REF.freq.y, A.ALT.freq.x,A.ALT.freq.y, B.ALT.freq.x, B.ALT.freq.y, C.ALT.freq.x,C.ALT.freq.y,D.ALT.freq.x,D.ALT.freq.y,A.B.suc.ALT.freq, A.B.suc.REF.freq, C.D.res.ALT.freq, C.D.res.REF.freq, Suc.AB.vs.Res.CD.REF.freq, A.vs.C.REF.freq, A.vs.D.REF.freq, B.vs.D.REF.freq, B.vs.C.REF.freq))

#reorder
d.v2 <- d.v2[order(-d.v2$FreqDiff_Susceptible_AB_vs_Resistant_CD),]


#write.csv(d.v2, file = "results/indels/GWSS_RNASeq1.snpEff.matrix_high.CombinedResults_R_indelfreqdifferences.csv")


#filter
d.v2.filt <- filter(d.v2, FreqDiff_Susceptible_AB_vs_Resistant_CD > 0.5 | FreqDiff_Tulare_Susceptible_A_vs_GeneralBeale_Resistant_C > 0.625 | FreqDiff_Tulare_Susceptible_A_vs_Tulare_Resistant_D > 0.625 | FreqDiff_Temecula_Susceptible_B_vs_GeneralBeale_Resistant_C > 0.625 | FreqDiff_Temecula_Susceptible_B_vs_Tulare_Resistant_D > 0.625)


#write.csv(d.v2.filt, file = "results/indels/GWSS_RNASeq1.snpEff.matrix_high.CombinedResults_R_indelFreqsDifferences.short.csv")
