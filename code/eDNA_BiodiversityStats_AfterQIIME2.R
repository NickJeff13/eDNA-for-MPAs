#load libraries
library("tidyverse")
library("ggplot2")

#set working directory to the dada2/qiime2 denoising output
setwd("~/eDNA/Musquash/Data/12S/dada2out_12S")

#load sample metadata and taxon list for ASVs
metadata <- read.delim("~/eDNA/Musquash/Data/12S/Musquash-12S-metadata_dada2.tsv", sep="\t") %>%  #metadata needs the first field to be "sampleid", all others can be custom
  rename(Site = sampleid) #rename the first column to Site
taxa <- read.delim("~/eDNA/Musquash/Data/12S/dada2out_12S/Musquash_12S_BlastClassifierTaxConsensus_ASV.tsv", sep="\t")  #manually split up the taxonomy into appropriate categories before uploading
readsPerSite <- read.delim("~/eDNA/Musquash/Data/12S/dada2out_12S/Musquash_12S_ASV_table_ReadsPerSite.tsv", sep="\t") #manually edit the ASV file to remove the first line of notes

#load reads per site and calculate total reads per ASV
readsPerSite <- readsPerSite %>% 
  left_join(select(taxa, ASV, Consensus, Taxon)) %>% #add in Taxon and Consensus score from taxa table
  relocate(Taxon, ASV, Consensus) %>% #Move Taxon, ASV, and Consensus columns to the beginning
  mutate_at(c(4:90), as.numeric) %>% #convert the read data to numeric
  rowwise() %>% #ensure the rest of the command is done across each row
  mutate(TotalReads = sum(c_across(4:90)), #sum all the reads per ASV across all samples
         OnePercTotalReads = round(TotalReads*0.01,0)) %>% #calculate 1% of total reads per ASV
  relocate(Taxon, ASV, Consensus, TotalReads, OnePercTotalReads)  #reorder the columns

#remove any ASV counts per site that are less than 1% total reads for the ASV
readsPerSite_clean <- readsPerSite %>% 
  rowwise() %>% #ensure the rest of the command is done across each row
  mutate_at(vars(c(6:92)), funs(ifelse(.<OnePercTotalReads,0,.))) #in cols 6:92, if the value is < OnePercTotalReads for that row, replace with 0

#Switch the table from wide to long
readsPerSite_clean_long <- readsPerSite_clean %>%
  gather(Site, Reads, X2A01S1:PCR_blank_002, factor_key=TRUE) %>% 
  relocate(Site,Taxon,Reads,ASV, Consensus, TotalReads, OnePercTotalReads) %>%
  filter(Reads !=0) %>%
  mutate(Site = str_replace(Site, "X","")) %>% 
  left_join(metadata)
write.csv(readsPerSite_clean_long, "Musquash_12S_ASV_table_ReadsPerSite_clean_long.csv")

###Make a table with the total number of filtered reads per site
TotalReadsPerSite <- readsPerSite_clean_long %>% 
  group_by(Site) %>% 
  summarize(TotalReads = sum(Reads))
write.csv(TotalReadsPerSite, "Musquash_12S_TotalReadsPerSite_filtered.csv")
