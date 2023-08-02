#load libraries
library("tidyverse")
library("ggplot2")
library("janitor")

#set working directory to the dada2/qiime2 denoising output
setwd("/mnt/sda2/eDNA/Musquash/Data/12S/dada2out_12S")

##before working with the reads table, make sure the taxonomy is finalized for the ASVs
BLASTOutput <- read.delim("/mnt/sda2/eDNA/Musquash/Data/12S/dada2out_12S/Musquash_12S_BlastClassifierTax_ASV.tsv", sep="\t",
                          col.names = c("ASV","Accession","PercIdent","AlignLength","Mismatches","GapOpens","QStart","QEnd","SStart","SEnd","Evalue","BitScore"))

BLASTConsensus <- read.delim("/mnt/sda2/eDNA/Musquash/Data/12S/dada2out_12S/Musquash_12S_BlastClassifierTaxConsensus_ASV.tsv", sep="\t")  #manually split up the taxonomy into separate columns categories before uploading

BLASTOutputTax <- BLASTOutput %>% 
  left_join(select(BLASTConsensus, ASV, Taxon)) %>% 
  rename(ConsensusTaxon = Taxon) %>% 
  relocate(ASV, ConsensusTaxon) 
write_delim(BLASTOutputTax, "/mnt/sda2/eDNA/Musquash/Data/12S/dada2out_12S/Musquash_12S_BlastClassifierConsensusTax_ASV_working.tsv", delim="\t")

##Manually edit the output to confirm taxonomy and BLAST unassigned ASVs, then load final list of ASV and Taxon
#This involves checking the PercIdent for each match, finding the lowest common taxon for matches >97%
BLASTOutput_FinalTax <- read.delim("/mnt/sda2/eDNA/Musquash/Data/12S/dada2out_12S/Musquash_12S_BlastClassifierConsensusTax_ASV_FINAL.tsv", sep="\t", col.names = c("ASV","Taxon"))

#load sample metadata ASVs
metadata <- read.delim("/mnt/sda2/eDNA/Musquash/Data/12S/Musquash-12S-metadata_dada2.tsv", sep="\t") %>%  #metadata needs the first field to be "sampleid", all others can be custom
  rename(Site = sampleid) #rename the first column to Site
readsPerSite_raw <- read.delim("/mnt/sda2/eDNA/Musquash/Data/12S/dada2out_12S/Musquash_12S_ASV_table_ReadsPerSite.tsv", sep="\t") #manually edit the ASV file to remove the first line of notes

#load reads per site and calculate total reads per ASV
readsPerSite <- readsPerSite_raw %>% 
  left_join(select(BLASTOutput_FinalTax, ASV, Taxon)) %>% #add in Taxon and Consensus score from taxa table
  relocate(Taxon, ASV) %>% #Move Taxon, ASV, and Consensus columns to the beginning
  mutate_at(c(4:89), as.numeric) %>% #convert the read data to numeric
  rowwise() %>% #ensure the rest of the command is done across each row
  mutate(TotalReadsASV = sum(c_across(4:89)), #sum all the reads per ASV across all samples
         OnePercTotalReadsASV = round(TotalReadsASV*0.01,0)) %>% #calculate 1% of total reads per ASV
  relocate(Taxon, ASV, TotalReadsASV, OnePercTotalReadsASV) #reorder the columns
  
###Make a table with the total number of filtered reads per site
TotalReadsPerSite <- as.data.frame(addmargins(readsPerSite[,5:91], c(1), sum)[4981,]) %>% 
  t() %>% 
  as.data.frame() %>% 
  rownames_to_column(.,"Site") %>% 
  mutate(Site = str_replace(Site, "X","")) %>% 
  rename(TotalReadsSite = V1) %>% 
  mutate(OnePercTotalReadsSite = round(TotalReadsSite*0.01,0))

#remove any ASV counts per site that are less than 1% total reads for the ASV
readsPerSite_ReadFiltered <- readsPerSite %>% 
  rowwise() %>% #ensure the rest of the command is done across each row
  mutate_at(vars(c(5:91)), funs(ifelse(.<OnePercTotalReadsASV,0,.))) #in cols 6:92, if the value is < OnePercTotalReads for that row, replace with 0

#Switch the table from wide to long1
readsPerSite_ReadFiltered_long <- readsPerSite_ReadFiltered %>%
  gather(Site, Reads, X2A01S1:PCR_blank_002, factor_key=TRUE) %>% 
  filter(Reads !=0) %>%
  mutate(Site = str_replace(Site, "X","")) %>% 
  left_join(metadata) %>% 
  left_join(TotalReadsPerSite) %>% 
  relocate(Site,Taxon,Reads,ASV,TotalReadsASV,OnePercTotalReadsASV, TotalReadsSite, OnePercTotalReadsSite)
write.csv(readsPerSite_ReadFiltered_long, "Musquash_12S_ASV_table_ReadsPerSite_ReadFiltered_long.csv")

#####remove any ASV counts per site that are less than 1% total reads for that SITE####
####ISSUE - losing a lot of species! many taxa are just being completely kicked out
#Skip this filtering for now
# readsPerSite_ReadSiteFiltered_long <- readsPerSite_ReadFiltered_long %>% 
#   group_by(Site,Taxon) %>% #group data by site and Taxon
#   mutate(ReadsPerTaxonPerSite = sum(Reads)) %>% #compute total reads per Site per Taxon
#   relocate(Site, Taxon,Reads, ASV, TotalReadsASV,OnePercTotalReadsASV, TotalReadsSite, OnePercTotalReadsSite, ReadsPerTaxonPerSite) %>% #reorder the columns
#   rowwise() %>% #ensure the rest of the command is done across each row
#   mutate_at(vars(Reads), funs(ifelse(ReadsPerTaxonPerSite<OnePercTotalReadsSite,0,Reads))) %>% #in Reads, if the value is < OnePercTotalReadsSite for that ASV, replace with 0
#   filter(Reads !=0) #remove any rows that have 0 reads
# write.csv(readsPerSite_ReadSiteFiltered_long, "Musquash_12S_TotalReadsPerSite_ReadSitefiltered_long.csv")
######

###Remove non-target taxa (humans, bacteria, etc)
AllTaxa <- unique(readsPerSite_ReadFiltered_long$Taxon) %>% #Create dataframe with all of the taxa
  as.data.frame() %>% 
  rename(Taxon = ".") 
# write.csv(AllTaxa, "Musquash_12S_ASV_AllTaxa.csv")
#add in the common names manually, then re-import
AllTaxa <- read.csv("Musquash_12S_ASV_AllTaxa.csv")

#filter readsPerSite_ReadFiltered_long to remove any of the Taxa to remove
readsPerSite_ReadFiltered_long_KeepTaxa <- readsPerSite_ReadFiltered_long %>% 
  filter(Taxon %in% subset(AllTaxa, AllTaxa$Keep == "Yes")$Taxon)

###Convert back to matrix of Site by ASV with reads, include Taxon
ReadsPerSite_Taxon_FinalFilter <- readsPerSite_ReadFiltered_long_KeepTaxa %>% 
  select(Site, Taxon, Reads) %>% #Keep only Site, Taxon, Reads
  pivot_wider(names_from = Site, values_from = Reads, values_fn = sum) %>% #add up the reads for any matching Taxon within a site
  replace(is.na(.),0)





