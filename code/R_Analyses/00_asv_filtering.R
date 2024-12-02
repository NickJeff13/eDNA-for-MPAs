library(dplyr)
library(tidyr)
library(tidyverse)

#simple function to filter reads in an ASV file - set to 0.01% now but can be changed to 0.1% etc
filter_low_read_species <- function(df) {
  numeric_df <-df[sapply(df,is.numeric)]
  total_reads <- sum(numeric_df, na.rm = TRUE)  # Total reads across species and sites
  species_sums <- rowSums(numeric_df, na.rm = TRUE)  # Sum of reads per species
  df_filtered <- df[species_sums >= 0.0001 * total_reads, ]  # Keep species with at least 1% of total reads
  return(df_filtered)
}

###########################2021 DATA###################################
##ESI

esi12s <- read.table(file = "data/2021Data/12s results/mergedspecies.tsv", header = T, sep="\t")
esi12s.filt <- esi12s %>%
  select(!ASV) %>%
  filter(pident>97 & Group %in% c("birds","bony fishes", "crustaceans", "gastropods", "hemichordates","isopods", "jellyfishes", "lancelets","ribbon worms","sea cucumbers","sea urchins","segmented worms","sharks & rays","starfish","whales & dolphins")) %>% 
  group_by(Species)

esi12s.filt.fish <- esi12s %>% select(!ASV) %>%
  filter(pident>98 & Group %in% "bony fishes") %>% 
  group_by(Species)

dim(esi12s.filt) #649 104
#remove Group, pident, and evalue columns next
esi12s.filt <- esi12s.filt[,-c(2:4)]
esi12s.filt.fish <- esi12s.filt.fish[,-c(2:4)]

#run this object in NMDS_plot.R script next

#16S Teleo
#ASV table
esi16s<-read.table(file = "data/2021Data/16s results/ESI16S_feature_table_export.tsv", header = T, sep = "\t")
#taxonomy
esi16s.taxa <-read.table(file = "data/2021Data/16s results/2116Sblast_results.tsv",header = F,sep="\t")
colnames(esi16s.taxa)<-c("ASV","NCBI","percentID", "evalue","length","species","group","score","commonname")
esi16s.asvs <-left_join(esi16s, esi16s.taxa, by="ASV")
dim(esi16s.asvs)
esi16s.filt<-esi16s.asvs %>% filter(percentID>98 & group=="bony fishes" & !species %in% c("Artediellus pacificus", "Liopsetta pinnifasciata", "Myzopsetta punctatissima")) #use this for NMDS and species accum plots


###############2022 Eastern Shore Data: 3 markers################################################################
#MiFishU 
esi22.12s <-read.table(file = "data/2022Data/ESI/MiFishU/ESI22_12S_feature_table_export.tsv", header = T, sep = "\t") %>% glimpse()
esi22.12s.taxa<-read.table("data/2022Data/ESI/MiFishU/12Sblast_results.tsv", sep="\t")

esi22.12s.merge <-left_join(esi22.12s, esi22.12s.taxa, by =c("ASV"="V1")) %>% distinct() %>% filter(V7 %in% c("bony fishes","whales & dolphins", "sharks & rays"), V3>97.99)

###############2022 SAB DATA#########################################################
sab16s<-read.table(data/2022Data/SAB/16S/SAB2216S_feature_table_FILTERED_forAPP.csv, sep=\t,header = T)
dim(sab16s) #2026 85

filt <- sab16s %>% filter(Confidence > 0.90) %>%
  group_by(Species) %>% 
  summarise(across(c(Sample.1,   Sample.10,  Sample.11,  Sample.12,  Sample.13,  Sample.14, 
              Sample.15,  Sample.16,  Sample.17, Sample.18 , Sample.19,  Sample.2,   Sample.20 , Sample.21,  Sample.22,
              Sample.23,  Sample.24,  Sample.25,  Sample.26,  Sample.27,  Sample.28,  Sample.29,  Sample.3  , Sample.30, 
              Sample.31,  Sample.32,  Sample.33,  Sample.34,  Sample.35 , Sample.36 , Sample.37,  Sample.38,  Sample.39, 
              Sample.4,   Sample.40,  Sample.41,  Sample.42,  Sample.43,  Sample.44 , Sample.45,  Sample.46,  Sample.47, 
              Sample.48,  Sample.49 , Sample.5,   Sample.50,  Sample.51,  Sample.52,  Sample.53,  Sample.54,  Sample.55, 
              Sample.56,  Sample.57,  Sample.58,  Sample.59,  Sample.6,   Sample.60,  Sample.61,  Sample.62,  Sample.63 ,
              Sample.64,  Sample.65,  Sample.66,  Sample.67,  Sample.68,  Sample.69,  Sample.7,   Sample.70,  Sample.71, 
              Sample.72,  Sample.73,  Sample.74,  Sample.75,  Sample.76,  Sample.77,  Sample.78,  Sample.79,  Sample.8,  
              Sample.80,  Sample.81, Sample.82, Sample.83,  Sample.9),sum)) %>%
  as.data.frame()

dim(filt) #1768 85
write.table(x = filt, file = "data/2022Data/SAB/16S/FilteredASVtable.txt",quote = F, sep = "\t")


#############COI##############
#Merging taxonomy and read count tables, and filtering taxonomy

#load in our two data tables - these paths will only work if you open the Courtney-Trask github R project, or set the working directory to your local github folder
tax <- read.table(file = "data/2022Data/SAB/COI/SABrdp5.output", header = F, sep = "\t") #here the \t is short for tab separation
tax<- rename(tax, OTU.ID= V1) #this is just so there is a matching column between our taxon and asvs files so we can merge them easier (OTU.ID)

asvs <- read.table(file="data/2022Data/SAB/COI/COI_feature_table_export.tsv", header = T, sep="\t")

taxatable <- merge(x= asvs, y= tax, by = "OTU.ID")  #this will reorganize the new table by OTU alphabetically

colnames(taxatable)
#Clean up this table and remove columns we don't need - here we are 'selecting' the columns from "taxatable" to keep, AND filtering by 0.90 probability of species being kept, AND removing rows with very few ASV counts

taxtable.final <- taxatable %>% select(OTU.ID, Sample.1,   Sample.10,  Sample.11,  Sample.12,  Sample.13,  Sample.14, 
                                       Sample.15,  Sample.16,  Sample.17, Sample.18 , Sample.19,  Sample.2,   Sample.20 , Sample.21,  Sample.22,
                                       Sample.23,  Sample.24,  Sample.25,  Sample.26,  Sample.27,  Sample.28,  Sample.29,  Sample.3  , Sample.30, 
                                       Sample.31,  Sample.32,  Sample.33,  Sample.34,  Sample.35 , Sample.36 , Sample.37,  Sample.38,  Sample.39, 
                                       Sample.4,   Sample.40,  Sample.41,  Sample.42,  Sample.43,  Sample.44 , Sample.45,  Sample.46,  Sample.47, 
                                       Sample.48,  Sample.49 , Sample.5,   Sample.50,  Sample.51,  Sample.52,  Sample.53,  Sample.54,  Sample.55, 
                                       Sample.56,  Sample.57,  Sample.58,  Sample.59,  Sample.6,   Sample.60,  Sample.61,  Sample.62,  Sample.63 ,
                                       Sample.64,  Sample.65,  Sample.66,  Sample.67,  Sample.68,  Sample.69,  Sample.7,   Sample.70,  Sample.71, 
                                       Sample.72,  Sample.73,  Sample.74,  Sample.75,  Sample.76,  Sample.77,  Sample.78,  Sample.79,  Sample.8,  
                                       Sample.80,  Sample.81,  Sample.82, Sample.83,  Sample.9, V12, V15, V27, V29) %>% filter(V12 %in% c("Nematoda","Arthropoda","Echinodermata","Porifera","Cnidaria","Mollusca","Annelida", "Platyhelminthes", "Nemertea","Chordata","Kinorhyncha","Brachiopoda","Ctenophora_combe_jellies","Bryozoa","Chaetognatha","Rotifera", "Hemichordata","Priapulida","Hemichordata", "Tardigrada") & V29 > 0.95) %>% as.data.frame()


coi.filt <- taxtable.final %>%
  group_by(V27) %>% 
  summarise(across(c(Sample.1,   Sample.10,  Sample.11,  Sample.12,  Sample.13,  Sample.14, 
                     Sample.15,  Sample.16,  Sample.17, Sample.18 , Sample.19,  Sample.2,   Sample.20 , Sample.21,  Sample.22,
                     Sample.23,  Sample.24,  Sample.25,  Sample.26,  Sample.27,  Sample.28,  Sample.29,  Sample.3  , Sample.30, 
                     Sample.31,  Sample.32,  Sample.33,  Sample.34,  Sample.35 , Sample.36 , Sample.37,  Sample.38,  Sample.39, 
                     Sample.4,   Sample.40,  Sample.41,  Sample.42,  Sample.43,  Sample.44 , Sample.45,  Sample.46,  Sample.47, 
                     Sample.48,  Sample.49 , Sample.5,   Sample.50,  Sample.51,  Sample.52,  Sample.53,  Sample.54,  Sample.55, 
                     Sample.56,  Sample.57,  Sample.58,  Sample.59,  Sample.6,   Sample.60,  Sample.61,  Sample.62,  Sample.63 ,
                     Sample.64,  Sample.65,  Sample.66,  Sample.67,  Sample.68,  Sample.69,  Sample.7,   Sample.70,  Sample.71, 
                     Sample.72,  Sample.73,  Sample.74,  Sample.75,  Sample.76,  Sample.77,  Sample.78,  Sample.79,  Sample.8,  
                     Sample.80,  Sample.81, Sample.82, Sample.83,  Sample.9),sum)) %>%
  as.data.frame()

dim(coi.filt) #124 84
write.table(x = coi.filt, file = "data/2022Data/SAB/COI/FilteredASVtable.txt",quote = F, sep = "\t")


##########################################################################################################################################################
################################################### Musquash 2023 data ###################################################################################
##########################################################################################################################################################

#Read in blast data and ASV table for 2023 12S data first
blasts <- read.csv("data/Musquash/2023/12Sblast_results_rarereads_filtered.csv",header = F,sep = "\t")
colnames(blasts)<-c("ASV","NCBIname","percentmatch","evalue","length","species","taxongroup","commonname")
head(blasts)

asvs <- read.csv("data/Musquash/2023/Musquash2023_filtered_table_biom/Musquash2023_feature_table_filtered.csv", header = T, sep="\t")
head(asvs)
dim(asvs)
dim(blasts)

#only use 'distinct()' if there are multiple blast hits per ASV
asv_taxa <- left_join(asvs, blasts, by="ASV") %>% distinct() %>%  filter(taxongroup %in% c("bony fishes","birds","carnivores","rodents","whales & dolphins","starfish","insectivores","sharks & rays"))
write.table(asv_taxa, file = "data/Musquash/2023/Musquash2023_TaxonTable_Filtered.tsv",sep = "\t",row.names = F)


###############################################################################################
#############################Beach Seining 2023################################################
###############################################################################################
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

m <-read.csv("data/2023Seining/12S/Seining_ASV_TaxonTable_Filtered_FishOnly.csv", header = T) %>% glimpse()
                                       
mm <- m %>% filter(PercentID > 98)

mmm <-filter_low_read_species(mm)


#Now the 2023 COI seining data
seining_coi <- read.table("data/2023Seining/COI-LERAYXT/ESI2023_featuretable_export.tsv", header = F)
inverts <- read.delim("data/2023Seining/COI-LERAYXT/rdp.output",header = F)

#merge the taxonomy from rdp.output and the feature table, and start by removing all mentions of bacteria
seine.inverts <-full_join(seining_coi, inverts, by = "V1") %>% filter(!grepl("bacteria", V6.y, ignore.case=TRUE))


###################################
#Once we have our filtered ASV tables, move to the NMDS scripts to make these plots