#GWSS SNP analysis between resistant & susceptible populations
#Adding gene annotations to hit list of SNPS
#Following GWSS_SNP_freq.R script
#Author: Cassie Ettinger
#email: cassandra.ettinger@ucr.edu

#load libraries
library(tidyverse)
library(magrittr)


#In terminal

#first get gene IDs from SNPs of interest
#cat GWSS_RNASeq1.snpEff.matrix_high.CombinedResults_R_AlleleFreqsDifferences.short.csv | cut -d ',' -f 7 > GWSS.CombinedResults_R_AlleleFreqsDifferences.short.GeneIDs.txt

#remove quotes
#sed 's/\"//g' GWSS.CombinedResults_R_AlleleFreqsDifferences.short.GeneIDs.txt > GWSS.CombinedResults_R_AlleleFreqsDifferences.short.GeneIDs.v2.txt

#replace orignal with quotes
#mv GWSS.CombinedResults_R_AlleleFreqsDifferences.short.GeneIDs.v2.txt GWSS.CombinedResults_R_AlleleFreqsDifferences.short.GeneIDs.txt

#mv to data folder

#make header for file
#head -n 1 Homalodisca_vitripennis_A6A7A9_masurca_v1.annotations.txt > HV_A6A7A9_masurca_v1.annotations.high.short.txt

# for  line in $(cat GWSS.CombinedResults_R_AlleleFreqsDifferences.short.GeneIDs.txt);
# do grep $line Homalodisca_vitripennis_A6A7A9_masurca_v1.annotations.txt | head -1 >> HV_A6A7A9_masurca_v1.annotations.high.short.txt;
# done

#reminder to gzip big annotation file to save space

#many duplicates - because multiple SNPs in same gene?
#cat HV_A6A7A9_masurca_v1.annotations.high.short.txt | sort | uniq -u > HV_A6A7A9_masurca_v1.annotations.high.short.sort.txt

#Now we have a reduced annotation file that we can easily match to our SNPS in R! YAY!

GWSS_annot <- read.delim("data/HV_A6A7A9_masurca_v1.annotations.high.short.txt")

SNP_targets <- read.csv("results/GWSS_RNASeq1.snpEff.matrix_high.CombinedResults_R_AlleleFreqsDifferences.short.csv")

#having trouble merging - keep getting duplicates 
#also issues with one gene name, fixing as follows
SNP_targets$GENE[SNP_targets$GENE == "LIG3"] <- "HOVITM_112457"

#remove some columns
GWSS_annot <- subset(GWSS_annot, select = -c(Notes, Alias.Synonyms, EC_number, InterPro, EggNog, COG, GO.Terms, gDNA, mRNA, CDS.transcript, Translation))

GWSS_SNP_annot <- inner_join(SNP_targets, GWSS_annot, by = c("GENE" = "GeneID", "CHROM" = "Contig"))


#using distinct to only get unique rows 
GWSS_SNP_annot.uniq <- GWSS_SNP_annot %>% distinct()


#write.csv(GWSS_SNP_annot.uniq, file = "results/GWSS_SNP_annotations.snpEff.high.freqdiff.short.csv")





### Add pfam gene names ###

# downloaded pfam_hits from funannoate folder
# opened in Excel to make into tab delim and unnecessary columns
pfam_names <- read.delim('data/pfam_hits.txt', header= FALSE)

#remove duplicate rows
pfam_names.uniq <- pfam_names %>% distinct()

#remove from R space
rm(pfam_names)

#assign header name to column for merge
pfam_names.uniq$GENE <- pfam_names.uniq$V3
pfam_names.uniq$Pfam.Ann <- pfam_names.uniq$V1
pfam_names.uniq$Pfam.ID <- pfam_names.uniq$V2

pfam_names.uniq <- subset(pfam_names.uniq, select = -c(V1, V2, V3))

#merge - will end up with duplicates bc some have multiple pfam hits
GWSS_SNP_annot.pfam <- left_join(GWSS_SNP_annot.uniq, pfam_names.uniq)

#remove from R space
rm(pfam_names.uniq)

#Group by SNP and collapse PFAM annotations (b/c there are multiple annotations per gene)
GWSS_SNP_annot.pfam.v2 <- GWSS_SNP_annot.pfam %>%  
  group_by(CHANGEDNA) %>% 
  summarise(PFAM.annot = str_c(Pfam.Ann, collapse = "; "))  %>%
  ungroup()

#merge new pfam annotations back to original SNP df 
GWSS_SNP_annot.pfam.annot <- left_join(GWSS_SNP_annot.uniq, GWSS_SNP_annot.pfam.v2)

#Remove fron R space
rm(GWSS_SNP_annot.pfam.v2)

#save output, note we are overwriting early output here so be careful! 
#write.csv(GWSS_SNP_annot.pfam.v2, file = "results/GWSS_SNP_annotations.snpEff.high.freqdiff.short.csv")
