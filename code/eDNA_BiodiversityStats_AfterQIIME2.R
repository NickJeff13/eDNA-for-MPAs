##########load libraries and RData#######
library("adegenet")
library("factoextra")
library("ggmap")
library("ggplot2")
library("ggrepel")
library("janitor")
library("paletteer")
library("raster")
library("RColorBrewer")
library("Rcpp")
library("rnaturalearth")
library("rnaturalearthdata")
library("rnaturalearthhires")
library("sf")
library("tidyverse")
library("vegan")
library("viridis")
library("wesanderson")

#set working directory to the dada2/qiime2 denoising output
setwd("C:/Users/VANWYNGAARDENMA/Documents/eDNA/Musquash/12S")
load("~/eDNA/Musquash/12S/Musquash_12S_workspace.RData")

#########Load and Process Data########
##before working with the reads table, make sure the taxonomy is finalized for the ASVs
BLASTOutput <- read.delim("Musquash_12S_BlastClassifierTax_ASV.tsv", sep="\t",
                          col.names = c("ASV","Accession","PercIdent","AlignLength","Mismatches","GapOpens","QStart","QEnd","SStart","SEnd","Evalue","BitScore"))

BLASTConsensus <- read.delim("Musquash_12S_BlastClassifierTaxConsensus_ASV.tsv", sep="\t")  #manually split up the taxonomy into separate columns categories before uploading

BLASTOutputTax <- BLASTOutput %>% 
  left_join(select(BLASTConsensus, ASV, Taxon)) %>% 
  rename(ConsensusTaxon = Taxon) %>% 
  relocate(ASV, ConsensusTaxon) 
# write_delim(BLASTOutputTax, "Musquash_12S_BlastClassifierConsensusTax_ASV_working.tsv", delim="\t")

##Manually edit the output to confirm taxonomy and BLAST unassigned ASVs, then load final list of ASV and Taxon
#This involves checking the PercIdent for each match, finding the lowest common taxon for matches >97%
BLASTOutput_FinalTax <- read.delim("Musquash_12S_BlastClassifierConsensusTax_ASV_FINAL.tsv", sep="\t", col.names = c("ASV","Taxon"))

#load sample metadata ASVs
metadata <- read.delim("Musquash-12S-metadata_dada2.tsv", sep="\t") %>%  #metadata needs the first field to be "sampleid", all others can be custom
  rename(Site = Sample) #rename the first column to Site
readsPerSite_raw <- read.delim("Musquash_12S_ASV_table_ReadsPerSite.csv", sep=",") #manually edit the ASV file to remove the first line of notes

#load reads per site and calculate total reads per ASV
readsPerSite <- readsPerSite_raw %>% 
  left_join(select(BLASTOutput_FinalTax, ASV, Taxon)) %>% #add in Taxon and Consensus score from taxa table
  relocate(Taxon, ASV) %>% #Move Taxon, ASV, and Consensus columns to the beginning
  mutate_at(c(3:89), as.numeric) %>% #convert the read data to numeric
  rowwise() %>% #ensure the rest of the command is done across each row
  mutate(TotalReadsASV = sum(c_across(3:89)), #sum all the reads per ASV across all samples
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
# write.csv(readsPerSite_ReadFiltered_long, "Musquash_12S_ASV_table_ReadsPerSite_ReadFiltered_long.csv")

###remove any ASV counts per site that are less than 1% total reads for that SITE##
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


###Remove non-target taxa (humans, bacteria, etc)
AllTaxa <- unique(readsPerSite_ReadFiltered_long$Taxon) %>% #Create dataframe with all of the taxa
  as.data.frame() %>% 
  rename(Taxon = ".") 
# write.csv(AllTaxa, "Musquash_12S_ASV_AllTaxa.csv")
#add in the common names manually and whether they should be kept (Common Name, Keep) then re-import
AllTaxa <- read.csv("Musquash_12S_ASV_AllTaxa.csv")

#filter readsPerSite_ReadFiltered_long to remove any of the Taxa to remove
readsPerSite_ReadFiltered_long_KeepTaxa <- readsPerSite_ReadFiltered_long %>% 
  filter(Taxon %in% subset(AllTaxa, AllTaxa$Keep == "Yes")$Taxon) 

#check how many ASVs each site has
ASVs_per_Site <- as.data.frame(table(readsPerSite_ReadFiltered_long_KeepTaxa$Site))

#remove the extraction blank - 1 ASV with only 1 read
readsPerSite_ReadFiltered_long_KeepTaxa <- readsPerSite_ReadFiltered_long_KeepTaxa %>% 
  filter(Site != "Extraction_blank_002")

#keep only bony fish in the long format for other analysis
readsPerSite_ReadFiltered_long_Fish <- readsPerSite_ReadFiltered_long_KeepTaxa %>% 
  filter(!Taxon %in% c("Phocidae","Phocoena phocoena","Squalus")) 

###Convert back to matrix of Site by ASV with reads, include Taxon
ReadsPerSite_Taxon_FinalFilter <- readsPerSite_ReadFiltered_long_KeepTaxa %>% 
  select(Site, Taxon, Reads) %>% #Keep only Site, Taxon, Reads
  pivot_wider(names_from = Taxon, values_from = Reads, values_fn = sum) %>% #add up the reads for any matching Taxon within a site
  replace(is.na(.),0) %>% 
  column_to_rownames(var = "Site")

##Keep only bony fish for further analysis
ReadsPerSite_Taxon_FinalFilter_Fish <- ReadsPerSite_Taxon_FinalFilter %>% 
  select(-c("Phocidae","Phocoena phocoena","Squalus"))

##Export Final Taxon per site (include common name) to make a table in Excel
ReadsPerSite_Taxon_FinalFilter_export <- as.data.frame(t(ReadsPerSite_Taxon_FinalFilter)) %>% 
  rownames_to_column("Taxon") %>% 
  # mutate_if(is.numeric, ~1 * (. > 0)) %>% #change all numbers greater than 0 to 1
  left_join(.,AllTaxa) %>% 
  relocate(Taxon,Common.Name) %>% 
  rename(CommonName=Common.Name)
write.csv(ReadsPerSite_Taxon_FinalFilter_export,"ReadsPerSite_Taxon_FinalFilter_12S_13Sept2023.csv" )



#######Site Map######
metadata <- read.delim("Musquash-12S-metadata_dada2.tsv", sep="\t") %>%  #metadata needs the first field to be "sampleid", all others can be custom
  rename(Site = Sample) #rename the first column to Site

#change lat/long to decimal degrees
metadata <- metadata %>% 
  arrange(site) %>% 
  mutate(Lat = as.numeric(str_replace(Lat, " N","")),
         Long = as.numeric(str_replace(Long, " W",""))) %>% 
  filter(!is.na(Lat)) %>% 
  mutate(Long = -Long) %>% 
  mutate(site = factor(site, levels = c("M01","M02","M04","M06","M07","M08","M11","M15","M16","M19",
                                           "M20","M22","M24","M26","M28","M30","M32",
                                           "2A01","2A10","2A12",
                                           "2B01","2B02","2B03",
                                           "301","302","303")))

metadata_plot <- metadata %>% 
  mutate(Label=site) %>% 
  distinct(site, .keep_all=T)

#Following Ryan's code from github
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
musquash <- read_sf("//ent.dfo-mpo.ca/ATLShares/Science/CESD/HES_MPAGroup/Data/Shapefiles/Musquash/Musquash_MPA_region.shp")%>%
  st_transform(latlong)
border <- read_sf("//ent.dfo-mpo.ca/ATLShares/Science/CESD/HES_MPAGroup/Data/Shapefiles/Coastline/Land_AtlCanada_ESeaboardUS.shp")%>%
  st_transform(latlong)

#colour palette for labels
pal <- c("2A"="magenta4","2B"="darkorange3","3"="firebrick3",Musq="black")

#Plot map
Musquash_SiteMap <- ggplot()+
  geom_sf(data=border)+
  geom_sf(data=musquash%>%
            filter(!is.na(ZONE)),
          aes(fill=ZONE),
          col="black")+
  scale_fill_viridis(discrete=T,direction = -1)+
  scale_colour_manual(values=pal,
                      guide="none") +
  geom_label(data = metadata_plot,
                   aes(x = Long,
                       y = Lat,
                       colour=zone),
             label = metadata_plot$Label,
             label.padding = unit(0.15, "lines"), # Rectangle size around label
             label.size = 0.35,
             nudge_y=ifelse(metadata_plot$Label=="2B03",0.0004,0))+
  labs(title="Musquash Marine Protected Area",fill="")+
  theme_bw()+
  theme(legend.position="right", 
        legend.title = element_blank(), 
        legend.text = element_text(size=28),
        legend.key = element_blank(),
        axis.text = element_text(size = 28), 
        axis.title = element_blank(), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin = unit(c(5,5,7,5), "mm"),
        panel.grid = element_blank())  +
  xlim(c(-66.34056,-66.22902))+#Long
  ylim(c(45.13733,45.19496))+#Lat
  annotate(geom = "text", x = -66.23, y = 45.139, size = 6, label = "Bay of\nFundy", colour="grey40")+
  annotate(geom = "text", x = -66.308, y = 45.182, size = 6, label = "Musquash\nRiver", colour="grey40")
ggsave(Musquash_SiteMap, 
       file = "Musquash_12S_SiteMap.pdf", 
       height = 10, 
       width = 16, 
       units = "in")

#######Taxonomy Summaries#########
ReadsPerSite_Taxon_FinalFilter_Fish_long <- ReadsPerSite_Taxon_FinalFilter_Fish %>% 
  mutate(Site = rownames(.)) %>% 
  pivot_longer(cols=c(1:ncol(ReadsPerSite_Taxon_FinalFilter_Fish)),names_to = "Species", values_to = "Reads") %>%    #add up the reads for any matching Taxon within a site
  filter(Reads != 0) %>% 
  left_join(metadata, by = "Site")

Taxon_by_Zone <- ReadsPerSite_Taxon_FinalFilter_Fish_long %>% 
  select(Species, zone) %>% 
  distinct(Species, .keep_all=T)




######Distance Matrices#######
#convert the reads matrix to a presence absence matrix
ReadsPerSite_FINAL_PresenceAbsence <- ReadsPerSite_Taxon_FinalFilter_Fish %>% 
  mutate_if(is.numeric, ~1 * (. > 0))#change all numbers greater than 0 to 1
ReadsPerSite_FINAL_PresenceAbsence <- na.omit(ReadsPerSite_FINAL_PresenceAbsence) #remove any rows with missing values


#Jaccard Distance
Jaccard <- vegdist(ReadsPerSite_FINAL_PresenceAbsence, 
                   method="jaccard", #jaccard is presence/absence, Bray-curtis is abundance
                   binary=TRUE) #presence absence standardization using decostand

BrayCurtis <- vegdist(ReadsPerSite_Taxon_FinalFilter_Fish,
                      method="bray",
                      binary=F)

########## Jaccard PCoA##########
###########generate PCoA using Jaccard distance
Musquash_12S_PCoA_Jaccard <- dudi.pco(Jaccard,
                             scannf = F,
                             nf=25)
summary(Musquash_12S_PCoA_Jaccard)
Musquash_12S_PCoA_Jaccard #view details

###Plot
#get dataframe ready for plot
Musquash_12S_PCoA_Jaccard_SiteData <- Musquash_12S_PCoA_Jaccard$li #site data
Musquash_12S_PCoA_Jaccard_SiteData$Site <- row.names(Musquash_12S_PCoA_Jaccard_SiteData)
Musquash_12S_PCoA_Jaccard_SiteData <- Musquash_12S_PCoA_Jaccard_SiteData %>% 
  left_join(metadata, by="Site") %>% 
  mutate(Sample = Site) %>% 
  mutate(zone = case_when(site %in% c("M01","M02","M04","M06","M07","M08","M11","M15","M16","M19") ~ "Zone 1",
                          site %in% c("M20","M22","M24","M26","M28","M30","M32","2A01","2A12") ~ "Zone 2A",
                          site %in% c("2B02","2B03") ~ "Zone 2B",
                          site %in% c("301","302","303") ~ "Zone 3")) %>% 
  relocate(Sample, type, location, site, zone,Lat,Long) %>%
  mutate(zone=factor(zone)) %>%
  arrange(zone)

####plot###
Musquash_12S_PCoA_Jaccard_baseplot <- ggplot(Musquash_12S_PCoA_Jaccard_SiteData, 
                                         aes(x=A1, 
                                             y=A2, 
                                             label=site))
Musquash_12S_PCoA_Jaccard_plot <- Musquash_12S_PCoA_Jaccard_baseplot +
  geom_hline(linewidth=1, colour = "black", yintercept = 0) + #horizontal line through the origin
  geom_vline(linewidth=1, colour = "black", xintercept = 0)  +#vertical line through the origin
  geom_point(aes(fill = site,
                 shape = site),
             colour="grey25",
             size=5, 
             alpha = 0.8,
             position=position_jitter(width=0, height=0.008)) + #transparency is alpha, individual points, coloured by Population
  xlab("\nPC1 (22.4%)")   +#label for the xaxis
  ylab("PC2 (15.3%)\n")   + #label for the yaxis
  scale_fill_manual(values=c(rep("#FDE725FF",times=5), #Zone1 FDE725FF
                               rep("#D2E21BFF",times=5), #Zone1
                               rep("#54C568FF",times=5), #Zone2A Musq 35B779FF
                               rep("#1F988BFF",times=4), #Zone2A
                               rep("#39568CFF",times=2), #Zone2B 31688EFF
                               rep("#440154FF",times=3))) + #Zone3 440154FF
  scale_shape_manual(values = c(rep(c(21,22,23,24,25),times=3),
                                21,22,23,24,
                                21,22,
                                21,22,23))+ #shapes for surface and bottom samples, circle 21, square22, triangle24, diamond23
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size=28),
        legend.key = element_blank(), 
        axis.text = element_text(size = 28), 
        axis.title=element_text(size = 35), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the legend, axes, etc - legend.key=element_blank() removes boxes from around the individual legend elements

ggsave(Musquash_12S_PCoA_Jaccard_plot, 
       file = "Musquash_12S_PCoA_biplot_Jaccard_13Sept2023.pdf", 
       height = 10, 
       width = 16, 
       units = "in")

########## Bray-Curtis PCoA##########
#######generate PCoA using Bray Curtis distance
Musquash_12S_PCoA_BC <- dudi.pco(BrayCurtis,
                                 scannf = F,
                                 nf=25)
summary(Musquash_12S_PCoA_BC)
Musquash_12S_PCoA_BC #view details

###Plot
#get dataframe ready for plot
Musquash_12S_PCoA_BC_SiteData <- Musquash_12S_PCoA_BC$li #site data
Musquash_12S_PCoA_BC_SiteData$Site <- row.names(Musquash_12S_PCoA_BC_SiteData)
Musquash_12S_PCoA_BC_SiteData <- Musquash_12S_PCoA_BC_SiteData %>% 
  left_join(metadata, by="Site") %>% 
  mutate(Sample = Site) %>% 
  mutate(zone = case_when(site %in% c("M01","M02","M04","M06","M07","M08","M11","M15","M16","M19") ~ "Zone 1",
                          site %in% c("M20","M22","M24","M26","M28","M30","M32","2A01","2A12") ~ "Zone 2A",
                          site %in% c("2B02","2B03") ~ "Zone 2B",
                          site %in% c("301","302","303") ~ "Zone 3")) %>% 
  relocate(Sample, type, location, site, zone,Lat,Long) %>%
  mutate(zone=factor(zone)) %>%
  arrange(zone)

####plot###
Musquash_12S_PCoA_BC_baseplot <- ggplot(Musquash_12S_PCoA_BC_SiteData, 
                                             aes(x=A1, 
                                                 y=A2, 
                                                 label=site))
Musquash_12S_PCoA_BC_plot <- Musquash_12S_PCoA_BC_baseplot +
  geom_hline(linewidth=1, colour = "black", yintercept = 0) + #horizontal line through the origin
  geom_vline(linewidth=1, colour = "black", xintercept = 0)  +#vertical line through the origin
  geom_point(aes(fill = site,
                 shape = site),
             colour="grey25",
             size=5, 
             alpha = 0.8,
             position=position_jitter(width=0, height=0.008)) + #transparency is alpha, individual points, coloured by Population
  xlab("\nPC1 (25.0%)")   +#label for the xaxis
  ylab("PC2 (16.8%)\n")   + #label for the yaxis
  scale_fill_manual(values=c(rep("#FDE725FF",times=5), #Zone1 FDE725FF
                             rep("#D2E21BFF",times=5), #Zone1
                             rep("#54C568FF",times=5), #Zone2A Musq 35B779FF
                             rep("#1F988BFF",times=4), #Zone2A
                             rep("#39568CFF",times=2), #Zone2B 31688EFF
                             rep("#440154FF",times=3))) + #Zone3 440154FF
  scale_shape_manual(values = c(rep(c(21,22,23,24,25),times=3),
                                21,22,23,24,
                                21,22,
                                21,22,23))+ #shapes for surface and bottom samples, circle 21, square22, triangle24, diamond23
  theme_bw() + #overall theme
  theme(legend.position = "right", 
        legend.title = element_blank(), 
        legend.text = element_text(size=28),
        legend.key = element_blank(), 
        axis.text = element_text(size = 28), 
        axis.title=element_text(size = 35), 
        panel.border = element_rect(linewidth =0.5),
        plot.margin=unit(c(5,5,7,5), "mm"),
        panel.grid.major=element_blank()) #changes the legend, axes, etc - legend.key=element_blank() removes boxes from around the individual legend elements

ggsave(Musquash_12S_PCoA_BC_plot, 
       file = "Musquash_12S_PCoA_biplot_BC_18Sept2023.pdf", 
       height = 10, 
       width = 16, 
       units = "in")



############Save Workspace#########
save.image("~/eDNA/Musquash/12S/Musquash_12S_workspace.RData")
