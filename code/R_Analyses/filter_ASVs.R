# Code to merge ASV tables and taxonomy from eDNA data processed in Qiime2, using dada2 and BLAST or RDP Classifier for taxonomy.

# Load libraries ----------------------------------------------------------

library(dplyr)
library(tidyr)
library(tidyverse)

#simple function to filter reads in an ASV file - set to 0.01% now but can be changed to 0.1% etc
filter_low_reads <- function(df) {
  numeric_df <-df[sapply(df,is.numeric)]
  total_reads <- sum(numeric_df, na.rm = TRUE)  # Total reads across species and sites
  species_sums <- rowSums(numeric_df, na.rm = TRUE)  # Sum of reads per species
  df_filtered <- df[species_sums >= 0.0001 * total_reads, ]  # Keep species with at least 0.01% of total reads
  return(df_filtered)
}

### Note: if exporting an ASV table from QIIME2 and using dada2, it will comment out the first two lines. You can read these in with header=F or edit the ASV tables before loading into R.

# 2021 Data ---------------------------------------------------------------

##ESI
### Update May 2025 - reprocessed 2021 data so created a 'NEW' folder for these results

esi12s<-read.table(file = "data/2021Data/NEW/12S/ESI2021_12S_feature_table_export.tsv", header = T, sep = "\t")
#taxonomy
esi12s.taxa <-read.table(file = "data/2021Data/NEW/12S/12Sblast_1results.tsv",header = F,sep="\t")
colnames(esi12s.taxa)<-c("ASV","NCBI","percentID", "evalue","length","species","group","commonname")

esi21.12s.merge <- left_join(esi12s, esi12s.taxa, by =c("OTU.ID"="ASV"))  %>% filter(group %in% c("bony fishes","whales & dolphins", "sharks & rays", "birds"), percentID > 97.99)



dim(esi21.12s.merge) #2464 114
#remove Group, pident, and evalue columns next

esi12s.filt.fish <- filter_low_reads(esi21.12s.merge)
esi12s.filt.fish <- esi12s.filt.fish %>% 
  dplyr::select(-c(NCBI,evalue,length,group,commonname, OTU.ID)) %>% 
  relocate(species) %>%
  relocate(percentID, .after=species)

write.csv(x = esi12s.filt.fish, file = "data/ESI2021_offshore_12S_filtered.csv", quote = F, row.names = F)
#run this object in NMDS_plot.R script next

#16S Teleo
#ASV table
esi16s<-read.table(file = "data/2021Data/NEW/16S/ESI2021_16S_feature_table_export.tsv", header = T, sep = "\t")
#taxonomy
esi16s.taxa <-read.table(file = "data/2021Data/NEW/16S/16Sblast_1results.tsv",header = F,sep="\t")
colnames(esi16s.taxa)<-c("ASV","NCBI","percentID", "evalue","length","species","group","commonname")

#can use inner_join or left_join for merging all ASV and taxonomy tables, as we will filter them by species/phylum next anyway. left_join will result in some rows having NA as not all ASVs get an assigned taxonomy, but inner_join will only keep rows in the ASV table that get a taxonomy as well 

esi21.16s.merge <-left_join(esi16s, esi16s.taxa, by=c("OTU.ID"="ASV"))
dim(esi21.16s.merge)
esi16s.filt <-esi21.16s.merge %>% filter(percentID>97.9 & group=="bony fishes" & !species %in% c("Artediellus pacificus", "Liopsetta pinnifasciata", "Myzopsetta punctatissima","Myoxocephalus polyacanthocephalus", "Gobio gobio", "Pholis laeta", "Psettichthys melanostictus","Platichthys environmental sample", "Hemitripterus villosus","Pseudopleuronectes yokohamae", "Sebastes steindachneri")) #use this for NMDS and species accum plots

esi16s.filt.fish <- filter_low_reads(esi16s.filt)
esi16s.filt.fish <- esi16s.filt.fish %>% 
  dplyr::select(-c(NCBI,evalue,length,group,commonname, OTU.ID)) %>% 
  relocate(species) %>%
  relocate(percentID, .after=species)

#write csv for Kayley for GOTeDNA
write.csv(x = esi16s.filt.fish, file = "data/ESI2021_offshore_16S_filtered.csv", quote = F, row.names = F)


# 2022 Eastern Shore Data - 3 markers -------------------------------------

## 12S MiFishU 
esi22.12s <-read.table(file = "data/2022Data/ESI/MiFishU/ESI22_12S_feature_table_export.tsv", header = T, sep = "\t") %>% glimpse()
esi22.12s.taxa <- read.csv("data/2022Data/ESI/MiFishU/12Sblast_results_updated.csv", header=F)
#uses only the first taxon from the blast results, which have 3 or 5 species listed
esi22.12s.taxa <- esi22.12s.taxa[!duplicated(esi22.12s.taxa$V1), ]


esi22.12s.merge <- left_join(esi22.12s, esi22.12s.taxa, by =c("ASV"="V1"))  %>%
  filter(V7 %in% c("bony fishes","whales & dolphins", "sharks & rays"), V3>97.9) %>%
  dplyr::select(-c(ASV, V2, V4,V5,V7, V8)) %>% 
  relocate(V6) %>%
  relocate(V3, .after=V6) %>%
  rename(Species = V6, PercentID=V3)

esi22.12s.merge$Species <- gsub("Ammodytes marinus", "Ammodytes sp.", esi22.12s.merge$Species)

write.csv(x = esi22.12s.merge, file = "data/2022Data/ESI/MiFishU/GOTeDNA_ESI2022_filteredASVs.csv", quote = F, row.names = F)


## 16S Fish
esi22.16s <- read.table("data/2022Data/ESI/16S/ESI22_16S_feature_table_export.tsv", header = T, sep="\t") %>% glimpse()
esi22.16s.taxa <- read.table("data/2022Data/ESI/16S/16Sblast_results.tsv", header=F, sep="\t")

esi22.16s.merge <- left_join(esi22.16s, esi22.16s.taxa, by =c("ASV"="V1"))  %>% 
  filter(V7 %in% c("bony fishes","whales & dolphins", "sharks & rays"), V3>97.9) %>%
  dplyr::select(-c(ASV, V2, V4,V5,V7, V8)) %>% 
  relocate(V6) %>%
  relocate(V3, .after=V6) %>%
  rename(Species = V6, PercentID=V3)

write.csv(x = esi22.16s.merge, file = "data/2022Data/ESI/16S/GOTeDNA_ESI2022_16SfilteredASVs.csv", quote = F, row.names = F)


## COI LerayXT primer
esi22.coi <- read.table("data/2022Data/ESI/LerayXT/ESI22_COI_feature_table_export.tsv", header = T, sep="\t")
esi22.coi.taxa <- read.table("data/2022Data/ESI/LerayXT/ESI2022.rdp.output", sep="\t") %>% glimpse()

esi22.coi.merge <- left_join(esi22.coi, esi22.coi.taxa, by=c("ASV"="V1")) %>% 
  filter(V29>0.97, V12 %in% c("Arthropoda","Platyhelminthes","Chordata","Annelida",
                              "Mollusca","Nematoda","Rhodophyta","Gastrotricha","Chlorophyta",
                              "Echinodermata","Brachiopoda","Porifera","Cnidaria","Nemertea",
                              "Haptophyta","Streptophyta","Hemichordata","Bryozoa","Ctenophora_comb_jellies",
                              "Tardigrada","Rotifera", "Chaetognatha","Prasinodermophyta")) %>% 
  select(!starts_with(c("ENEG","EXT","PCRB"))) %>%
  rename(Phylum=V12, Class=V15, Species=V27)

esi22.coi.filt <- esi22.coi.merge %>% 
  dplyr::select(-c(ASV,V2:V26,V28)) %>% 
  relocate(Species) %>%
  relocate(V29, .after=Species) %>%
  rename(Confidence=V29)

write.csv(x = esi22.coi.filt, file = "data/2022Data/ESI/LerayXT/GOTeDNA_ESI2022_COIfilteredASVs.csv", quote = F, row.names = F)


# 2022 SAB Data -----------------------------------------------------------

## 12S MiFishU

sab2022.12s.asvs <-read.csv(file = "data/2022Data/SAB/MiFishU/SAB2022_12S_ASVtable.csv", header = T) %>% glimpse()
#taxonomy
sab22.12s.taxa <-read.csv(file = "data/2022Data/SAB/MiFishU/12S_blastresults.csv",header = F)

sab22.12s.merge <- left_join(sab2022.12s.asvs, sab22.12s.taxa, by =c("ASV"="V1"))  %>% 
  filter(V7 %in% c("bony fishes","whales & dolphins", "sharks & rays"), V3 > 97.9)


dim(sab22.12s.merge) #884 99
#remove Group, pident, and evalue columns next

sab22.12s.filt <- filter_low_reads(sab22.12s.merge)
sab22.12s.filt <- sab22.12s.filt %>% 
  rename(species=V6, percentID=V3) %>%
  dplyr::select(-c(ASV,V2, V4, V5, V7,V8)) %>% 
  relocate(species) %>%
  relocate(percentID, .after=species)

write.csv(x = sab22.12s.filt, file = "data/2022Data/SAB/MiFishU/GOTeDNA_SAB2022_12S_filtered.csv", quote = F, row.names = F)
#run this object in NMDS_plot.R script next

## 16S fish
sab22.16s<-read.table("data/2022Data/SAB/16S/SAB2216S_feature_table_FILTERED_forAPP.csv", sep="\t", header = T)
dim(sab22.16s) #2026 85

#this groups species but for GOTeDNA and MCT app we probably want non-aggregated (i.e., ASV level)
# filt <- sab16s %>% filter(Confidence > 0.90) %>%
#   group_by(Species) %>% 
#   summarise(across(c(Sample.1,   Sample.10,  Sample.11,  Sample.12,  Sample.13,  Sample.14, 
#               Sample.15,  Sample.16,  Sample.17, Sample.18 , Sample.19,  Sample.2,   Sample.20 , Sample.21,  Sample.22,
#               Sample.23,  Sample.24,  Sample.25,  Sample.26,  Sample.27,  Sample.28,  Sample.29,  Sample.3  , Sample.30, 
#               Sample.31,  Sample.32,  Sample.33,  Sample.34,  Sample.35 , Sample.36 , Sample.37,  Sample.38,  Sample.39, 
#               Sample.4,   Sample.40,  Sample.41,  Sample.42,  Sample.43,  Sample.44 , Sample.45,  Sample.46,  Sample.47, 
#               Sample.48,  Sample.49 , Sample.5,   Sample.50,  Sample.51,  Sample.52,  Sample.53,  Sample.54,  Sample.55, 
#               Sample.56,  Sample.57,  Sample.58,  Sample.59,  Sample.6,   Sample.60,  Sample.61,  Sample.62,  Sample.63 ,
#               Sample.64,  Sample.65,  Sample.66,  Sample.67,  Sample.68,  Sample.69,  Sample.7,   Sample.70,  Sample.71, 
#               Sample.72,  Sample.73,  Sample.74,  Sample.75,  Sample.76,  Sample.77,  Sample.78,  Sample.79,  Sample.8,  
#               Sample.80,  Sample.81, Sample.82, Sample.83,  Sample.9),sum)) %>%
#   as.data.frame()
sab22.16s.filt <- sab22.16s %>% select(-OTU.ID) %>%
  mutate(Species=na_if(Species, "")) %>%
  filter(!is.na(Species), Confidence>0.97, !Species==" Iluocoetes fimbriatus")

dim(sab22.16s.filt) #637 84


write.csv(x = sab22.16s.filt, file = "data/2022Data/SAB/16S/GOTeDNA_SAB202216SFilteredASVtable.csv",quote = F, row.names = F)


#COI general diversity 
#Merging taxonomy and read count tables, and filtering taxonomy

#load in our two data tables - these paths will only work if you open the Courtney-Trask github R project, or set the working directory to your local github folder
sab22.coi.taxa <- read.table(file = "data/2022Data/SAB/COI/SABrdp5.output", header = F, sep = "\t") #here the \t is short for tab separation
sab22.coi.taxa<- rename(sab22.coi.taxa, OTU.ID= V1) #this is just so there is a matching column between our taxon and asvs files so we can merge them easier (OTU.ID)

sab22.coi.asvs <- read.table(file="data/2022Data/SAB/COI/SAB22_COI_feature_table_export.tsv", header = T, sep="\t")

sab22.coi.merge <- left_join(sab22.coi.asvs, sab22.coi.taxa, by ="OTU.ID")  %>%
  filter(V29>0.97, V12 %in% c("Arthropoda","Platyhelminthes","Chordata","Annelida","Mollusca","Nematoda","Rhodophyta","Gastrotricha","Chlorophyta","Echinodermata","Brachiopoda","Porifera","Cnidaria","Nemertea","Haptophyta","Streptophyta","Hemichordata","Bryozoa","Ctenophora_comb_jellies","Tardigrada","Rotifera", "Chaetognatha","Prasinodermophyta")) %>% select(!starts_with(c("ENEG","EXT","PCRB"))) %>%
  rename(Phylum=V12, Class=V15, Species=V27)


#Clean up this table and remove columns we don't need - here we are 'selecting' the columns from "taxatable" to keep, AND filtering by 0.90 probability of species being kept, AND removing rows with very few ASV counts

sab22.coi.merge.filt <- sab22.coi.merge %>% filter(Phylum %in% c("Nematoda","Arthropoda","Echinodermata","Porifera","Cnidaria","Mollusca","Annelida", 
                                                              "Platyhelminthes", "Nemertea","Chordata",
                                                              "Kinorhyncha","Brachiopoda","Ctenophora_combe_jellies","Bryozoa","Chaetognatha",
                                                              "Rotifera", "Hemichordata","Priapulida","Hemichordata", "Tardigrada") & V29 > 0.95) %>% 
  select(-c(OTU.ID, V2:V11, Phylum, V13:V14, Class, V16:V26, V28)) %>%
  relocate(Species) %>%
  relocate(V29, .after=Species) %>%
  as.data.frame()


# coi.filt <- taxtable.final %>%
#   group_by(V27) %>% 
#   summarise(across(c(Sample.1,   Sample.10,  Sample.11,  Sample.12,  Sample.13,  Sample.14, 
#                      Sample.15,  Sample.16,  Sample.17, Sample.18 , Sample.19,  Sample.2,   Sample.20 , Sample.21,  Sample.22,
#                      Sample.23,  Sample.24,  Sample.25,  Sample.26,  Sample.27,  Sample.28,  Sample.29,  Sample.3  , Sample.30, 
#                      Sample.31,  Sample.32,  Sample.33,  Sample.34,  Sample.35 , Sample.36 , Sample.37,  Sample.38,  Sample.39, 
#                      Sample.4,   Sample.40,  Sample.41,  Sample.42,  Sample.43,  Sample.44 , Sample.45,  Sample.46,  Sample.47, 
#                      Sample.48,  Sample.49 , Sample.5,   Sample.50,  Sample.51,  Sample.52,  Sample.53,  Sample.54,  Sample.55, 
#                      Sample.56,  Sample.57,  Sample.58,  Sample.59,  Sample.6,   Sample.60,  Sample.61,  Sample.62,  Sample.63 ,
#                      Sample.64,  Sample.65,  Sample.66,  Sample.67,  Sample.68,  Sample.69,  Sample.7,   Sample.70,  Sample.71, 
#                      Sample.72,  Sample.73,  Sample.74,  Sample.75,  Sample.76,  Sample.77,  Sample.78,  Sample.79,  Sample.8,  
#                      Sample.80,  Sample.81, Sample.82, Sample.83,  Sample.9),sum)) %>%
#   as.data.frame()

dim(sab22.coi.merge.filt) #313  85
write.csv(x = sab22.coi.merge.filt, file = "data/2022Data/SAB/COI/GOTeDNA_SAB2022_COI.filtered.csv",quote = F, row.names = F)


# 2023 ESI Perley Data ----------------------------------------------------
#12S Fish
esi23.12s.perl <-read.table("data/2023Perley/ESI/MiFishU/ESI23_12S_feature_table_export.tsv", header = T, sep="\t") %>% glimpse()

esi23.12s.perl.taxa <-read.table("data/2023Perley/ESI/MiFishU/12Sblast_results.tsv", header = F, sep="\t") %>% glimpse()

esi23.12s.perl.merge <- left_join(esi23.12s.perl, esi23.12s.perl.taxa, by=c("ASV"="V1"))

esi23.12s.per.merge2 <- esi23.12s.perl.merge %>% 
  filter(V3>97.9 & V7 %in% c("bony fishes","whales & dolphins","sharks & rays"))

esi23.12s.perl.merge3 <- esi23.12s.per.merge2 %>%
  select(-c(ASV, V2, V4, V5, V7, V8)) %>%
  rename(species=V6, percentID=V3) %>%
  relocate(species) %>%
  relocate(percentID, .after=species)

write.csv(x = esi23.12s.perl.merge3, file="data/2023Perley/ESI/MiFishU/GOTeDNA_Perley2023_12S_formatted.csv", quote = F, row.names = F)

#COI Invertebrates
esi23.coi.perl <- read.table("data/2023Perley/ESI/COI/ESIPer23_COI_feature_table_export.tsv", header = T, sep = "\t") %>% glimpse()

esi23.coi.perl.rdp <- read.table("data/2023Perley/ESI/COI/ESIPerl.2023.rdp.output", header = F, sep="\t") %>% glimpse()

esi23.coi.perl.merge <- left_join(esi23.coi.perl, esi23.coi.perl.rdp, by=c("ASV"="V1"))

esi23.coi.perl.filt <- esi23.coi.perl.merge %>% 
  filter(V29>0.97, V12 %in% c("Arthropoda","Platyhelminthes","Chordata","Annelida",
                              "Mollusca","Nematoda","Rhodophyta","Gastrotricha","Chlorophyta",
                              "Echinodermata","Brachiopoda","Porifera","Cnidaria",
                              "Nemertea","Haptophyta","Hemichordata","Bryozoa",
                              "Ctenophora_comb_jellies","Tardigrada","Rotifera", 
                              "Chaetognatha","Kinorhyncha","Acanthocephala_thorny-headed_worms"),!V15=="Insecta") %>% 
                                          select(-c(ASV, V2:V26, V28)) %>%
                                          as.data.frame()

esi23.coi.perl.filt2 <- filter_low_reads(esi23.coi.perl.filt) %>% 
  rename(species=V27, confidence=V29) %>%
  relocate(species) %>%
  relocate(confidence, .after=species)
  
dim(esi23.coi.perl.filt2) #227 58

write.csv(x = esi23.coi.perl.filt2, file="data/2023Perley/ESI/COI/GOTeDNA_Perley2023_COI.csv", quote=F, row.names = F)

# 2023 SAB Perley Data ----------------------------------------------------

#12S Fish
sab23.12s.perl <-read.table("data/2023Perley/SAB/MiFishU/SAB2023_12S_feature_table_export.tsv", header = T, sep="\t") %>% glimpse()

sab23.12s.perl.taxa <-read.table("data/2023Perley/SAB/MiFishU/12Sblast_results.tsv", header = F, sep="\t") %>% glimpse()

sab23.12s.perl.merge <- left_join(sab23.12s.perl, sab23.12s.perl.taxa, by=c("OTU.ID"="V1"))

sab23.12s.per.merge2 <- sab23.12s.perl.merge %>% 
  filter(V3>97.9 & V7 %in% c("bony fishes","whales & dolphins","sharks & rays", "chimaeras")) %>%
  as.data.frame()

sab23.12s.perl.merge3 <- filter_low_reads(sab23.12s.per.merge2) %>%
  select(-c(OTU.ID, V2, V4, V5, V7, V8)) %>%
  rename(species=V6, percentID=V3) %>%
  relocate(species) %>%
  relocate(percentID, .after=species)

write.csv(x = sab23.12s.perl.merge3, file="data/2023Perley/SAB/MiFishU/GOTeDNA_SABPerley2023_12S_formatted.csv", quote=F, row.names = F)

#COI Invertebrates
sab23.coi.perl <- read.table("data/2023Perley/SAB/Leray/SAB2023_COI_feature_table_export.tsv", header = T, sep = "\t") %>% glimpse()

sab23.coi.perl.rdp <- read.table("data/2023Perley/SAB/Leray/rdp.output", header = F, sep="\t") %>% glimpse()

sab23.coi.perl.merge <- left_join(sab23.coi.perl, sab23.coi.perl.rdp, by=c("OTU.ID"="V1"))

sab23.coi.perl.filt <- sab23.coi.perl.merge %>% 
  select(-c(OTU.ID, V2:V11,V13,V22, V25, V28)) %>%
           filter(V26>0.92, V12 %in% c("Arthropoda","Platyhelminthes","Chordata","Annelida","Mollusca","Nematoda",
                                       "Rhodophyta","Gastrotricha","Chlorophyta","Echinodermata","Brachiopoda","Porifera",
                                       "Cnidaria","Nemertea","Haptophyta","Hemichordata","Bryozoa","Ctenophora_comb_jellies",
                                       "Tardigrada","Rotifera", "Chaetognatha","Kinorhyncha","Acanthocephala_thorny-headed_worms"),
                  !V15=="Insecta") %>%
           rename(Phylum=V12, Class=V15, species=V27) %>% as.data.frame()

sab23.coi.perl.filt2 <- filter_low_reads(sab23.coi.perl.filt) %>% 
  rename(confidence=V29) %>%
  relocate(species) %>%
  relocate(confidence, .after=species)

dim(sab23.coi.perl.filt2) #215 140

write.csv(x = sab23.coi.perl.filt2, file="data/2023Perley/SAB/Leray/GOTeDNA_SABPerley2023_COI.csv", quote=F, row.names = F)

# 2023 Musquash Data ------------------------------------------------------

#Read in blast data and ASV table for 2023 12S data first
blasts <- read.csv("data/Musquash/2023/12Sblast_results_rarereads_filtered.csv", header = F, sep = "\t")
colnames(blasts)<-c("ASV","NCBIname","percentmatch","evalue","length","species","taxongroup","commonname")
head(blasts)

asvs <- read.csv("data/Musquash/2023/Musquash2023_filtered_table_biom/Musquash2023_feature_table_filtered.csv", header = T, sep="\t")
head(asvs)
dim(asvs)
dim(blasts)

#only use 'distinct()' if there are multiple blast hits per ASV
asv_taxa <- left_join(asvs, blasts, by="ASV") %>% distinct() %>%  filter(taxongroup %in% c("bony fishes","birds","carnivores","rodents","whales & dolphins","starfish","insectivores","sharks & rays"))
write.table(asv_taxa, file = "data/Musquash/2023/Musquash2023_TaxonTable_Filtered.tsv",sep = "\t",row.names = F)



# 2023 ESI Beach Seining --------------------------------------------------

#load in the fish 12S ASV table and blast taxonomy 

seining_asvs <- read.table("data/2023Seining/12S/Seining12S_feature_table_export.tsv",header = T) %>% glimpse()

fishies <- read.delim("data/2023Seining/12S/seining23_blast.tsv", header = F, sep="\t") %>% distinct(V1,.keep_all = T) #For fish we will just keep the top hit for each ASV, but some rules need to be followed 

#Rules for Fish taxonomy
#1. All Pholis are likely to be P. gunnelus. However, we'll just go with Pholis sp. for analyses
#2. Shorthorn and Grubby sculpin have very similar sequences are difficult to tell apart. Any Myoxocephalus should just be Myoxocephalus sp.
#3. Pacific and Atlantic herring are almost indistinguishable with 12S. Change all Pacific herring to Atlantic (Clupea harengus)

seining_spp <- left_join(seining_asvs, fishies, by=c("ID"="V1"))
seining_spp <- seining_spp  %>% dplyr::select(-contains("PCRBlank"),-V2,-V3,-V6,-V8,-V9,-V10, -V11, -V12, -V13,-V14, -V15, -V16,-V17, -V18, -V20,-V22)
#Write the output 
write.csv(seining_spp, file="data/2023Seining/12S/Seining_ASV_TaxonTable_Filtered.csv", row.names = F)

#did some processing by hand, just removing non-fish. Do some more filtering for low reads and low percent matches 

esi23.seine.12s.merge <-read.csv("data/2023Seining/12S/Seining_ASV_TaxonTable_Filtered_FishOnly.csv", header = T) %>% glimpse()
                                       
esi23.seine.12s.merge.filt <- esi23.seine.12s.merge %>% filter(PercentID > 97.99)

esi23.seine.12s.merge.filt2 <-filter_low_reads(esi23.seine.12s.merge.filt) %>% 
  select(-c(ID, Common, Evalue))

#Write csv for Kayley for GOTeDNA
write.csv(x = esi23.seine.12s.merge.filt2, file = "data/2023Seining/12S/GOTeDNA_ESI2023_Coastal_12S.csv", quote = F, row.names = F)

#Now the 2023 COI seining data
esi23.coi.asvs <-read.table("data/2023Seining/COI-LERAYXT/ESI2023_featuretable_export.tsv", header = T, sep="\t") %>% glimpse()
esi23.coi.taxa <-read.table("data/2023Seining/COI-LERAYXT/rdp.output", header = F, sep="\t") %>% select(c("V1","V12","V15", "V24","V26","V27","V29"))

esi23.coi.merge <-left_join(esi23.coi.asvs,esi23.coi.taxa, by=c("OTU.ID"="V1"))

#filter
esi23.coi.filt <- esi23.coi.merge %>% filter(V26>0.94, V12 %in% c("Arthropoda","Platyhelminthes","Chordata","Annelida","Mollusca","Nematoda","Rhodophyta","Gastrotricha","Chlorophyta","Echinodermata","Brachiopoda","Porifera","Cnidaria","Nemertea","Haptophyta","Hemichordata","Bryozoa","Ctenophora_comb_jellies","Tardigrada","Rotifera", "Chaetognatha")) %>% select(!starts_with(c("ENEG","EXT","PCRB"))) %>%
  rename(Phylum=V12, Class=V15, Species=V27, Confidence=V29)

esi23.coi.filt2 <- filter_low_reads(esi23.coi.filt) %>% filter(!Species %in% c("Sus_scrofa","Homo_sapiens")) %>%
  select(-c(OTU.ID,Phylum,Class,V24,V26)) %>%
  relocate(Species) %>%
  relocate(Confidence, .after=Species)

#Write csv for Kayley for GOTeDNA
write.csv(x = esi23.coi.filt2, file = "data/2023Seining/COI-LERAYXT/GOTeDNA_ESI2023_Coastal_COI.csv", quote = F, row.names = F)



# 2024 Coastal Seining Data -----------------------------------------------

## 12S
  esi24.12s.coast <-read.table("data/2024Seining/MiFishU/ESI12S_feature_table_export.tsv", header = T, sep="\t")   %>% glimpse()

  esi24.12s.coast.taxa <-read.table("data/2024Seining/MiFishU/12Sblast_results.tsv", header = F, sep="\t") %>%   glimpse()

  esi24.12s.coast.merge <- left_join(esi24.12s.coast, esi24.12s.coast.taxa, by=c("OTU.ID"="V1"))

  esi24.12s.coast.filt <- filter_low_reads(esi24.12s.coast.merge %>% 
    drop_na() %>%
    filter(V3>97.99 & V7 %in% c("bony fishes","whales & dolphins","sharks & rays")) %>%
    select(-c(OTU.ID,V2, V4, V5, V8, V7))) %>%
    relocate(V6) %>%
    relocate(V3, .after=V6)
  
  esi24.12s.coast.filt$V6 <- gsub("Clupea pallasii","Clupea harengus", esi24.12s.coast.filt$V6)
  esi24.12s.coast.filt$V6 <- gsub("Pholis crassispina","Pholis gunnellus", esi24.12s.coast.filt$V6)
  

write.csv(x = esi24.12s.coast.filt,"data/2024Seining/MiFishU/ESI2024_MiFish_TaxonomyMerged.csv", quote=F)

## COI
esi24.coi.coast <- read.table("data/2024Seining/COI/ESI2024_COI_feature_table_export.tsv", header = T, sep = "\t") %>% glimpse()

esi24.coi.coast.rdp <- read.table("data/2024Seining/COI/rdp.output", header = F, sep="\t") %>% glimpse()

esi24.coi.coast.merge <- left_join(esi24.coi.coast, esi24.coi.coast.rdp, by=c("ASV"="V1"))

#normally I filter our insects but keeping them for now for the inland Keji site we sampled 

esi24.coi.coast.filt <- filter_low_reads(esi24.coi.coast.merge %>% 
                                          select(-c(V2:V11,V13,V22, V25, V28, Undetermined,ASV)) %>%
                                          filter(V26>0.92, V12 %in% c("Arthropoda","Platyhelminthes","Chordata","Annelida","Mollusca","Nematoda","Rhodophyta",
                                                                      "Gastrotricha","Chlorophyta","Echinodermata","Brachiopoda","Porifera","Cnidaria",
                                                                      "Nemertea","Haptophyta","Hemichordata","Bryozoa","Ctenophora_comb_jellies","Tardigrada",
                                                                      "Rotifera", "Chaetognatha","Kinorhyncha","Acanthocephala_thorny-headed_worms")) %>%
                                          rename(Phylum=V12, Class=V15, Species=V27) %>% 
                                           as.data.frame())

  write.csv(x = esi24.coi.coast.filt, "data/2024Seining/COI/ESI2024_COI_TaxonomyMerged.csv", quote = F)
  
  esi24.coi.coast.filt2 <- esi24.coi.coast.filt %>%
                           select(-c(Phylum, V14, Class, V16, V17, V18, V19,
                                     V20, V21, V23, V24, V26)) %>%
    relocate(Species) %>%
    relocate(V29, .after = Species)
  
  write.csv(esi24.coi.coast.filt2, file= "data/2024Seining/COI/GOTeDNA_ESI2024_COI_formatted.csv", quote=F, row.names = F)
  

# St Anns Bank 2024 -------------------------------------------------------

## These also include some AZMP samples from Browns Bank, Northeast Channel, and the Gully 
# Run by ABL in winter 2024, 248_F + MiFishU-R for 12S and Leray for COI 
  
## 12S
  sab24.12s.asv <-read.table("data/2024Perley/12S/denoised/SAB2024_12S_filtered_table_biom/SAB2024_12S_feature_table_export.tsv", header = T, sep="\t")   %>% glimpse()
  
  sab24.12s.taxa <-read.table("data/2024Perley/12S/denoised/12Sblast_1results.tsv", header = F, sep="\t") %>%   glimpse()
  
  sab24.12s.merge <- left_join(sab24.12s.asv, sab24.12s.taxa, by=c("OTU.ID"="V1"))
  
  sab24.12s.filt <- filter_low_reads(sab24.12s.merge %>% 
                                             drop_na() %>%
                                             filter(V3>97.99 & V7 %in% c("bony fishes","whales & dolphins","sharks & rays")) %>%
                                             select(-c(OTU.ID,V2, V4, V5, V8, V7))) %>%
    relocate(V6) %>%
    relocate(V3, .after=V6)
  
  write.csv(x = sab24.12s.filt,"data/2024Perley/12S/SAB2024_Perley_12S_filtered.csv", row.names = F, quote=F)
  

## COI 
  sab24.coi.asv <- read.table("data/2024Perley/COI/SAB2024_COI_feature_table_export.tsv", header = T, sep = "\t") %>% glimpse()
  
  sab24.coi.rdp <- read.table("data/2024Perley/COI/rdp.output", header = F, sep="\t") %>% glimpse()
  
  sab24.coi.merge <- left_join(sab24.coi.asv, sab24.coi.rdp, by=c("OTU.ID"="V1"))
  
  sab24.coi.filt <- filter_low_reads(sab24.coi.merge %>% 
                                             select(-c(V2:V11,V13,V22, V25, V28, OTU.ID)) %>%
                                             filter(V29>0.979, V12 %in% c("Arthropoda","Platyhelminthes","Chordata","Annelida","Mollusca","Nematoda","Rhodophyta",
                                                                         "Gastrotricha","Chlorophyta","Echinodermata","Brachiopoda","Porifera","Cnidaria",
                                                                         "Nemertea","Haptophyta","Hemichordata","Bryozoa","Ctenophora_comb_jellies","Tardigrada",
                                                                         "Rotifera", "Chaetognatha","Kinorhyncha","Acanthocephala_thorny-headed_worms")) %>%
                                             rename(Phylum=V12, Class=V15, Species=V27) %>% 
                                             as.data.frame())
  
  write.csv(x = sab24.coi.filt, "data/2024Perley/COI/SAB2024_COI_TaxonomyMerged.csv", quote = F)
  
  sab24.coi.filt2 <- sab24.coi.filt %>%
    select(-c(Phylum, V14, Class, V16, V17, V18, V19,
              V20, V21, V23, V24, V26)) %>%
    relocate(Species) %>%
    relocate(V29, .after = Species)
  
  write.csv(sab24.coi.filt2, file= "data/2024Perley/COI/GOTeDNA_SAB2024_COI_formatted.csv", quote=F, row.names = F)
  
  
  
# Notes -------------------------------------------------------------------

#Once we have our filtered ASV tables, move to the NMDS or diversity scripts to make some plots
#Below is a list of invertebrate distributions for those found in these eDNA datasets but I am not sure they are accurate

#Copepods
#1. Pseudocalanus mimus = Pacific, change to Pseuodcalanus sp.1
#2. Eurytemora herdmani = Atlantic, keep as is
#3. Temora longicornis = Atlantic, keep as is
#4. Pseudocalanus acuspes = Arctic, Baltic, some North Pacific, change to Pseudocalanus sp. 2
#5. Calanus glacialis = everywhere in the north, keep as is
#6. Calanus hyperboreus = Arctic and north Atlantic, keep as is
#7. Pseudocalanus newmani = Primarily Pacific and ARctic, possible records in Atlantic but change to Pseudocalanus sp.3

#Non-copepod crustaceans
#1. Thysanoessa_inermis, krill = NW Atlantic, Gulf species, keep


#Sessile Inverts
#1. Boltenia echinata = Atlantic and Arctic, keep as is
#2. Hymeraphia_stellifera, a sponge = European, change to Hymeraphia sp. 
#3. Haliclona oculata, a sponge = North Atlantic, keep as is
#4. Ectyonopsis_pluridentata, sponge = South African species only, remove?
#5. Crella_elegans sponge = Mediterranean species only, remove?
#6. Obelia geniculata, a hydroid= found here, leave as is


#Annelids
#1. Pholoe minuta, a polychaete = Atlantic, keep as is
#2. Galathownia oculata, a polychaete = Atlantic, including N. Am and Europe
#3. Pherusa plumoisa = North and South Atlantic, keep as is
#4. Gattyana cirrhosa = North Atlantic and Arctic, keep as is
#5. Dodecaceria_concharum = European records, some NW Atlantic but possible species mis-IDs? change to Dodecaceria sp?
#6. Prionospio_steenstrupi = NW atlantic, keep as is
#7. Micronephthys_minuta = NW Atlantic, keep as is
#8. Laonice_cirrata = NW Atlantic, keep
#9. Pseudopotamilla_reniformis = North Atlantic, keep as is 
#10. Pectinaria_hyperborea = North Atlantic, but accepted genus is Cistenides. Change genus name
#11. Trochochaeta multisetosa = found here, leave as is


#Molluscs
#1. Macoma calcarea = North Atlantic, keep as is
#2. Eubranchus_olivaceus, a sea slug = North Pacific, rename to Eubranchus sp. 

#Diatoms, algae, and phytoplankton
#1. Micromonas pusilla, Chlorophyta = Atlantic, mainly European but some N. Am records
#2. Dasysiphonia_japonica = invasive species, known in our region
#3. Ahnfeltia_plicata = red algae, keep as is
#4. Euthora_cristata red algae = keep as is

