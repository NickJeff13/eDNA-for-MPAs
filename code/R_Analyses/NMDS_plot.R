# Code for running species accumulation curves, NMDS, and other diversity metrics on a whole bunch of eDNA metabarcoding data including species richness plots, NMDS plots, and species accumulation plots. Written by Nick.Jeffery@dfo-mpo.gc.ca 

# Load libraries ----------------------------------------------------------

library(vegan)
library(dplyr)
library(janitor)
library(ggplot2)
library(tidyverse)
library(ggalluvial) #make alluvial plots
library(patchwork) #stick plots together
library(eulerr) #for Venn diagrams
library(pals)

# NMDS theme for all plots ------------------------------------------------

nmdstheme <- theme_bw()+
  theme(#axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text=element_text(size=18))


# Load RData object -------------------------------------------------------
# Can load this object if saved at the bottom of this script
load("data/eDNA_NMDS_and_Diversity.RData")


# 2021 ESI Perley Data ----------------------------------------------------

#### metadata
esi21meta <-read.table("data/2021Data/metadata/2021-sample-metadata_ESIonly.tsv",header = T, sep = "\t")
esi21meta$sample.id<-gsub("-",".",esi21meta$sample.id)

#### 12S
glimpse(esi12s.filt.fish)
#esi12smat <- esi12s.filt %>% group_by(Species) %>% summarise(across(everything(), sum)) %>% data.frame()
esi12tt <- t(esi12s.filt.fish[,2:length(colnames(esi12s.filt.fish))])
colnames(esi12tt)<-esi12s.filt.fish$Species
esi12ttt<-esi12tt[rowSums(esi12tt[])>0,]

##### Make a barplot of taxa
tt<-pivot_longer(esi12s.filt.fish, cols=starts_with("Sample"))
  
  p1 <- ggplot()+geom_bar(data=tt%>%filter(value>2000), aes(x=Species, y=log(value)),stat="identity")+
    xlab(label = "")+
    ylab(label="12S Log(Read Count)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust=1), text=element_text(size=14))+
    ggtitle('2021');p1

#16S
glimpse(esi16s.filt)
spec.mat<-as.data.frame(t(esi16s.filt[,2:101]))
colnames(spec.mat)<-esi16s.filt$species
spec.mat<-spec.mat[rowSums(spec.mat)>0,]

#Make a barplot of taxa
tt<-pivot_longer(esi16s.filt, cols=starts_with("Sample"))

  p2 <- ggplot()+geom_bar(data=tt%>%filter(value>3000), aes(x=species, y=log(value)),stat="identity")+
    xlab(label = "Species")+
    ylab(label="16S Log(Read Count)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust=1), text=element_text(size=14));p2

#plot with patchwork
p1/p2

ggsave(filename = "ESI2021_FishBarplots_Combined.png",plot = last_plot(), device = "png", path = "figures/2021Results/", width = 10, height=8, dpi = 300, bg = "white")

# Run NMDS on the ungrouped 12S and 16S datasets - better to leave ungrouped as the ASVs contribute to beta diversity better this way 
esi16s.nmds.jac<-metaMDS(comm = spec.mat, distance = "jaccard", k=3, trymax=100)
esi16s.nmds.bray<-metaMDS(comm = spec.mat, distance = "bray", k=3, trymax=100)
esi12s.nmds.jac<-metaMDS(comm = esi12ttt, distance = "jaccard", k=3, trymax=100)
esi12s.nmds.bray<-metaMDS(comm = esi12ttt, distance = "bray", k=3, trymax=100)

#extract nmds scores for ggplot - here just swap the various NMDS objects to make data.scores
data.scores = as.data.frame(scores(esi12s.nmds.jac)$sites)
data.scores$Sample <- rownames(data.scores)
data.scores.metadat <-left_join(data.scores,esi21meta, by=c("Sample"="sample.id"))


species.scores <- as.data.frame(scores(esi12s.nmds.jac, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.data.per21 <- data.scores.metadat %>%
  as.data.frame() %>%
  group_by(surface) %>%
  slice(chull(x=NMDS1,y=NMDS2)) #for this dataset, there's no differentiation between bottom and surface


#plot it up
  p3 <- ggplot(data.scores.metadat %>% filter(surface !="BLANK"), aes(x = NMDS1, y = NMDS2)) + 
   geom_polygon(data=hull.data.per21 %>% filter(surface !="BLANK"),aes(x=NMDS1,y=NMDS2,fill=surface),color="black",alpha=0.30) +
    scale_fill_manual(values=c("#c43b3b", "cornflowerblue"))+
    #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
    geom_jitter(size = 4, aes(fill = surface, shape = surface), 
                width = 10, height = 10) +
    #scale_fill_viridis_c(option="viridis") + # Continuous color scale for depth
    scale_shape_manual(values = c(21, 24)) +         # Discrete shapes (e.g., circle and triangle)
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "none", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) + 
    labs(x = "NMDS1", y = "NMDS2")  + 
   geom_text(aes(x=Inf, y=Inf, vjust=30,hjust=1.1,label=paste("Stress =",round(esi12s.nmds.bray$stress,3),"k =",esi12s.nmds.bray$ndim)));p3

#extract nmds scores for ggplot - here just swap the various NMDS objects to make data.scores
data.scores = as.data.frame(scores(esi16s.nmds.jac)$sites)
data.scores$Sample <- rownames(data.scores)
data.scores.metadat <-left_join(data.scores,esi21meta, by=c("Sample"="sample.id"))


species.scores <- as.data.frame(scores(esi16s.nmds.jac, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.data.per21 <- data.scores.metadat %>%
  as.data.frame() %>%
  group_by(surface) %>%
  slice(chull(x=NMDS1,y=NMDS2)) #for this dataset, there's no differentiation between bottom and surface


  p4 <- ggplot(data.scores.metadat %>% filter(surface !="BLANK"), aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data=hull.data.per21 %>% filter(surface !="BLANK"),aes(x=NMDS1,y=NMDS2,fill=surface),color="black",alpha=0.30) +
    scale_fill_manual(values=c("#c43b3b", "cornflowerblue"))+
    #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
    geom_jitter(size = 4, aes(fill = surface, shape = surface), 
              width = 10, height = 10) +
    #scale_fill_viridis_c(option="viridis") + # Continuous color scale for depth
    scale_shape_manual(values = c(21, 24)) +         # Discrete shapes (e.g., circle and triangle)
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
    labs(x = "NMDS1", y = "NMDS2")+
    geom_text(aes(x=Inf, y=Inf, vjust=30,hjust=1.1,label=paste("Stress =",round(esi16s.nmds.bray$stress,3),"k =",esi16s.nmds.bray$ndim)));p4

p3 + p4

ggsave(filename = "ESI21_12S_and_16S_NMDS1_2_combined.png",plot = last_plot(), device = "png", path = "figures/2021Results/", width = 10, height=8, dpi = 320)

# Shannon and other diversity metrics for each sample and station

shan21.12s <- diversity(esi12ttt, index="shannon") %>% 
  as.data.frame()
colnames(shan21.12s)<-c("Sample","ShannonDiv")

shan21.16s <- diversity(spec.mat, index="shannon") %>%
  as.data.frame()
simps21.12s <- diversity(esi12ttt, index="simpson") %>% 
  as.data.frame()
simps21.16s <- diversity(spec.mat, index="simpson") %>% 
  as.data.frame()

# 2022 ESI Perley Data ----------------------------------------------------

# load metadata
esi22.meta <-read.csv("data/2022Data/2022-sample-metadata_ESI.csv", header = T, sep="\t") %>% glimpse()
esi22.meta <- esi22.meta[-1,] %>% 
  select(-c("temperature","ph","dissolvedoxy","salinity"))#remove the first row which is useless from QIIME2

esi22.meta$sample.id <- gsub("-",".",esi22.meta$sample.id) # replace dashes with . to match the ASV table

###### 12S Data ##
head(esi22.12s.merge)
esi22.12s.merge$V6 <- gsub("Clupea pallasii", "Clupea harengus", esi22.12s.merge$V6)
esi22.12s.merge$V6 <- gsub("Sebastes baramenuke", "Sebastes sp.", esi22.12s.merge$V6)
esi22.12s.merge$V6 <- gsub("Sebastes viviparus", "Sebastes sp.", esi22.12s.merge$V6)
esi22.12s.merge$V6 <- gsub("Ammodytes hexapterus", "Ammodytes sp.", esi22.12s.merge$V6)
esi22.12s.merge$V6 <- gsub("Ammodytes personatus", "Ammodytes sp.", esi22.12s.merge$V6)
esi22.12s.merge$V6 <- gsub("Pollachius virens", "Pollachius pollachius", esi22.12s.merge$V6)
esi22.12s.merge$V6 <- gsub("Pholis ornata", "Pholis gunnellus", esi22.12s.merge$V6)


#Make a barplot of taxa
tt<-pivot_longer(esi22.12s.merge, cols=starts_with("Sample"))


  p5 <- ggplot()+geom_bar(data=tt%>%filter(value>2000), aes(x=V6, y=log(value)),stat="identity")+
    xlab(label = "Species")+
    ylab(label="12S Log(Read Count)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust=1), text=element_text(size=14))+
    ggtitle('2022');p5


###### 16S Data ##
head(esi22.16s.merge)
esi22.16s.merge$V6 <- gsub("Sebastes mentella", "Sebastes sp.", esi22.16s.merge$V6)
esi22.16s.merge$V6 <- gsub("Gadus macrocephalus", "Gadus morhua", esi22.16s.merge$V6)
esi22.16s.merge$V6 <- gsub("Pholis laeta", "Pholis gunnellus", esi22.16s.merge$V6)
esi22.16s.merge$V6 <- gsub("Platichthys environmental sample", "Platichthys flesus", esi22.16s.merge$V6)
esi22.16s.merge$V6 <- gsub("Pollachius virens", "Pollachius pollachius", esi22.16s.merge$V6)
esi22.16s.merge$V6 <- gsub("Pleuronectes platessa", "Pleuronectinae", esi22.16s.merge$V6)


tt<-pivot_longer(esi22.16s.merge, cols=starts_with("Sample"))

  p6 <- ggplot()+geom_bar(data=tt%>%filter(value>2000, !V6 %in% c("Platichthys flesus","Myzopsetta punctatissima"  )), aes(x=V6, y=log(value)),stat="identity")+
   xlab(label = "Species")+
    ylab(label="16S Log(Read Count)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust=1), text=element_text(size=14));p6

#plot with patchwork
p5/p6

ggsave(filename = "ESI2022_FishBarplots_Combined.png",plot = last_plot(), device = "png", path = "figures/2022Results//", width = 10, height=8, dpi = 300, bg = "white")

#Run the NMDS on the 12S and 16S data
#12S

mat22.12s<-as.data.frame(t(esi22.12s.merge %>% select(contains("Sample."))))
colnames(mat22.12s)<-esi22.12s.merge$V6
mat22.12s<-mat22.12s[rowSums(mat22.12s)>0,]

#Remove blanks (Samples 6, 12, 13, 28, 53, 88-91) and one outlier sample (Sample 35)
esi22.12s2 <-mat22.12s[-c(4,5,17,23,59,60,62),]

#16S
glimpse(esi22.16s.merge)
mat22.16s<-as.data.frame(t(esi22.16s.merge %>% select(starts_with("Sample"))))
colnames(mat22.16s)<-esi22.16s.merge$V6
mat22.16s<-mat22.16s[rowSums(mat22.16s)>0,]

esi22.16s2<-mat22.16s[-c(4,13,29),]
# Run NMDS on the ungrouped 12S and 16S datasets - better to leave ungrouped as the ASVs contribute to beta diversity better this way 
esi16s.nmds.jac<-metaMDS(comm = esi22.16s2, distance = "jaccard", k=8, trymax=100)
esi16s.nmds.bray<-metaMDS(comm = esi22.16s2, distance = "bray", k=8, trymax=100)
esi12s.nmds.jac<-metaMDS(comm = esi22.12s2, distance = "jaccard", k=8, trymax=100)
esi12s.nmds.bray<-metaMDS(comm = esi22.12s2, distance = "bray", k=8, trymax=100)

#extract nmds scores for ggplot - here just swap the various NMDS objects to make data.scores
data.scores = as.data.frame(scores(esi12s.nmds.jac)$sites)
data.scores$Sample <- rownames(data.scores)
data.scores.metadat <-left_join(data.scores,esi22.meta, by=c("Sample"="sample.id"))


species.scores <- as.data.frame(scores(esi12s.nmds.jac, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.data.per22 <- data.scores.metadat %>%
  as.data.frame() %>%
  group_by(surface) %>%
  slice(chull(x=NMDS1,y=NMDS2)) #for this dataset, there's no differentiation between bottom and surface


#plot it up
  p7 <- ggplot(data.scores.metadat, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data=hull.data.per22,aes(x=NMDS1,y=NMDS2,fill=surface),color="black",alpha=0.30) +
    scale_fill_manual(values=c("#c43b3b", "cornflowerblue"))+
    #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
    geom_jitter(data=data.scores.metadat, size = 4, 
           aes(fill = surface, shape = surface))+
                                              #scale_fill_viridis_c(option="viridis") + # Continuous color scale for depth
    scale_shape_manual(values = c(21, 24)) +         # Discrete shapes (e.g., circle and triangle)
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "none", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
    labs(x = "NMDS1", y = "NMDS2")  + 
    geom_text(aes(x=Inf, y=Inf, vjust=30,hjust=1.1,label=paste("Stress =",round(esi12s.nmds.jac$stress,3),"k =",esi12s.nmds.jac$ndim)));p7

#extract nmds scores for ggplot - here just swap the various NMDS objects to make data.scores
data.scores = as.data.frame(scores(esi16s.nmds.jac)$sites)
data.scores$Sample <- rownames(data.scores)
data.scores.metadat <-left_join(data.scores,esi22.meta, by=c("Sample"="sample.id"))


species.scores <- as.data.frame(scores(esi16s.nmds.jac, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.data.per22 <- data.scores.metadat %>%
  as.data.frame() %>%
  group_by(surface) %>%
  slice(chull(x=NMDS1,y=NMDS2)) #for this dataset, there's no differentiation between bottom and surface


#plot it up - 16S
  p8 <- ggplot(data.scores.metadat, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data=hull.data.per22,aes(x=NMDS1,y=NMDS2,fill=surface),color="black",alpha=0.30) +
    scale_fill_manual(values=c("#c43b3b", "cornflowerblue"))+
    #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
    geom_jitter(data=data.scores.metadat, size = 4, 
              aes(fill = surface, shape = surface))+
    #scale_fill_viridis_c(option="viridis") + # Continuous color scale for depth
    scale_shape_manual(values = c(21, 24)) +         # Discrete shapes (e.g., circle and triangle)
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) + 
    labs(x = "NMDS1", y = "")  + 
    geom_text(aes(x=Inf, y=Inf, vjust=30,hjust=1.1,label=paste("Stress =",round(esi16s.nmds.jac$stress,3),"k =",esi16s.nmds.jac$ndim)));p8


p7+p8+plot_annotation(title="2022",tag_levels = "A",
                    theme = theme(plot.title = element_text(size = 20)))

ggsave(filename = "ESI2022_12s_16s_NMDS_Jaccard_Combined.png", plot = last_plot(), device = "png", path = "figures/2022Results/", width = 10, height=8, dpi = 300, bg = "white")

###### LerayXT COI Data ##
head(esi22.coi.merge)

#Make a barplot of taxa
ii <-pivot_longer(esi22.coi.merge, cols=starts_with("Sample")) %>% filter(!Species=="Homo_sapiens")
ii$Species <- gsub("Nothria_conchylega_CMC02","Nothria_conchylega", ii$Species)
ii$Species <- gsub("Bipalponephtys_neotena","Micronephthys_neotena", ii$Species)
ii$Species <- gsub("Euclymene_sp.","Euclymene_zonalis", ii$Species)

  p9 <- ggplot()+geom_bar(data=ii%>%filter(value>400), aes(x=Species, y=log(value)),stat="identity")+
    xlab(label = "Species")+
    ylab(label="COI Log(Read Count)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust=1), text=element_text(size=14));p9

  p9+plot_annotation(title = "2022",
                   theme = theme(plot.title = element_text(size = 20)))

ggsave(filename = "ESI2022_COI_Barplot.png",plot = last_plot(), device = "png", path = "figures/2022Results/", width = 10, height=8, dpi = 300, bg = "white")

#COI NMDS
spec.mat.coi<-as.data.frame(t(esi22.coi.merge %>% select(starts_with("Sample"))))
colnames(spec.mat.coi)<-esi22.coi.merge$Species
spec.mat.coi<-spec.mat.coi[rowSums(spec.mat.coi)>0,]

#Remove blanks - samples 6,12,13,28, 53, 88-91 for this dataset
esi22.coi.mm2 <-spec.mat.coi[-c(4,5,21,84),]
nmds.esi22.coi <- metaMDS(esi22.coi.mm2, distance = "jaccard", k=18 , trymax = 200, maxit=500)

#extract nmds scores for ggplot
data.scores = as.data.frame(scores(nmds.esi22.coi)$sites)
data.scores$Sample <- rownames(data.scores)

#merge with metadata by sample id column
data.scores.merge <- inner_join(data.scores,esi22.meta, by=c("Sample"="sample.id"))


species.scores <- as.data.frame(scores(nmds.esi22.coi, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.data.per22.coi <- data.scores.merge %>%
  as.data.frame() %>%
  group_by(surface) %>%
  slice(chull(x=NMDS1,y=NMDS2)) #for this dataset, there's no differentiation between bottom and surface


  p10 <- ggplot() + 
    geom_polygon(data=hull.data.per22.coi,aes(x=NMDS1, y=NMDS2, fill=surface),color="black",alpha=0.30) +
    scale_fill_manual(values=c("firebrick","cornflowerblue"),name="Water sample")+
    ggnewscale::new_scale_fill()+
    geom_point(data=data.scores.merge, aes(x = NMDS1, y = NMDS2, fill=as.numeric(depth),shape=surface),size = 4, colour="black")+ 
    scale_shape_manual(values=c(21,24),name="Water sample")+
    scale_fill_viridis_c(option="viridis",name="Depth (m)")+
    #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
          #legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "none", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          #legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
    labs(x = "NMDS1", y = "NMDS2")  + 
    geom_text(aes(x=Inf, y=Inf, vjust=55,hjust=1.2,label=paste("Stress =",round(nmds.esi22.coi$stress,3),"k =",nmds.esi22.coi$ndim)))+
    nmdstheme;p10

  p10 + plot_annotation(title = "Perley 2022 COI",
                                    theme = theme(plot.title = element_text(size = 20)))

ggsave(filename = "ESI22_Perley_COI_NMDS1_2_byDepth.png",plot = last_plot(), device = "png", path = "figures/2022Results/", width = 12, height=8, units = "in", dpi = 400, bg = "white")  



# 2022 SAB Perley Data ----------------------------------------------------


#read in our data table for species accummulation curves and NMDS plots
taxtable <- read.table("data/2022Data/SAB/COI/COI_FilteredASVtable.txt", header = T)

#Make a barplot of taxa
tt<-pivot_longer(taxtable, cols=starts_with("X"))

  p11 <- ggplot()+geom_bar(data=tt, aes(x=V27, y=log(value)),stat="identity")+
    xlab(label = "Species")+
    ylab(label="Log(Read Count)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.2), text=element_text(size=14));p11

ggsave(filename = "SAB2022_COI_barplot.png",plot = last_plot(), device = "png", path = "figures/", width = 10, height=8, units = "in", dpi = 400, bg = "white")

meta <- read.table("data/2022Data/SAB/R_2022-sample-metadata_SABonly.tsv", header = T)
rownames(meta)<-meta$watersample
meta<-meta[1:79,]
head(meta)


commat <- taxtable #%>% select(V26, Sample.1,  Sample.2, Sample.3,  Sample.5, Sample.6, Sample.7, Sample.8, Sample.9,Sample.10, Sample.12, Sample.13, Sample.14, Sample.16, Sample.17, Sample.18, Sample.20, Sample.21, Sample.22, Sample.24, Sample.25, Sample.26) #removed samples 4, 11, 15, 19, 23, 27 as these are the negatives

#Transpose the table for vegan
commat2<-t(commat[,2:length(colnames(commat))])
commat2<-commat2[rowSums(commat2[])>0,]
colnames(commat2)<-commat[,1]

groups <- read.csv("data/2022Data/SAB/SAB_COI_depthdata.csv", header = F)



sab.coi.nmds <-metaMDS(commat2, distance="bray", k=6, trymax = 100, maxit=500)
plot(sab.coi.nmds) #this is not very informative without labels!



ordiplot(sab.coi.nmds, type='n')
ordihull(sab.coi.nmds,groups=groups$V2,draw="polygon",col="cyan",label=F)
#orditorp(sab.coi.nmds,display="species",col="red",air=0.01)
orditorp(sab.coi.nmds,display="sites",   air=0.01,cex=0.75)

#extract nmds scores for ggplot
data.scores = as.data.frame(scores(sab.coi.nmds)$sites)
data.scores$Sample <- rownames(data.scores)
data.scores$Depth <- groups$V2
data.scores$Surface <- groups$V3


species.scores <- as.data.frame(scores(sab.coi.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 


  p12 <- ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
    geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
    geom_point(size = 4, aes(shape = Surface, colour = Depth))+ 
      theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) + 
    labs(x = "NMDS1", colour = "Depth", y = "NMDS2", shape = "Surface")  + 
    geom_text(aes(x=Inf, y=Inf, vjust=65,hjust=1.2,label=paste("Stress =",round(sab.coi.nmds$stress,3),"k =",sab.coi.nmds$ndim)))+  
    scale_colour_continuous(trans="reverse");p12


ggsave(filename = "SAB2022_COI_NMDS.png",plot = p12, device = "png", path = "figures/", width = 10, height=8, units = "in", dpi = 400, bg = "white")

#Accumulation curves
sp_list.coi <- taxtable$V27
raremin<-min(rowSums(commat2))

yy<-specaccum(commat2,method="exact", permutations = 1000)

tidy_specaccum <- function(x) {
  data.frame(
    site = x$sites,
    richness = x$richness,
    sd = x$sd)
}
yyy <- tidy_specaccum(yy)

  p13 <- ggplot() +
    geom_line(data=yyy, aes(x=site, y=richness), linewidth=2, color="firebrick") +
    geom_linerange(data=yyy,aes(x = site, ymin = richness - 2*sd, ymax = richness + 2*sd)) +
    ylim(0, NA)+
    ylab(label = "Species Richness")+
    xlab(label="Site")+
    theme_bw()+
    theme(text = element_text(size=20));p13

  ggsave(filename = "2022SAB_COI_Specaccum.png", plot = p13, device = "png", 
         path = "figures/", width = 10, height=8, units="in",dpi = 400, bg="white")

###Now do the same thing for 12S! 
######################################
fishdat <- read.table("data/2022Data/SAB/16S/16S_FilteredASVtable.txt", sep="\t",header = T)
s.groups <- read.csv("data/2022Data/SAB/SAB16S_GroupData.csv")
head(fishdat)
#get rid of row 1 which is non-fish animals
fishdat <- fishdat[-c(1,16),]

mm<-pivot_longer(fishdat, cols=starts_with("X"))


  p14 <- ggplot()+
    geom_bar(data=mm, aes(x=Species, y=log(value)),stat="identity")+
    xlab(label = "Species")+
    ylab(label="log(Read Count)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), text=element_text(size=18));p14

ggsave(filename = "2022SAB_16S_barplot.png", plot = last_plot(), device = "png", path = "figures/", width = 10, height=8, units="in",dpi = 400, bg="white")

fishdat.t<-t(fishdat[,2:length(colnames(fishdat))])
fishdat.t<-fishdat.t[rowSums(fishdat.t[])>0,]
colnames(fishdat.t)<-fishdat[,1]

sab.16s.nmds <-metaMDS(fishdat.t, distance="bray", k=6, trymax = 200, maxit=200)
plot(sab.16s.nmds)



data.scores.16s = as.data.frame(scores(sab.16s.nmds)$sites)

data.scores.16s$Sample <- rownames(data.scores.16s)
data.scores.16s$Depth <- s.groups$depth
data.scores.16s$Surface <- s.groups$surface

species.scores <- as.data.frame(scores(sab.16s.nmds, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 


ff = ggplot(data.scores.16s, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(shape = Surface, colour = Depth))+ 
  geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species))+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", colour = "Depth", y = "NMDS2", shape = "Surface")  + 
  scale_colour_continuous(trans="reverse") 

ff

ggsave(filename = "SAB2022_16S_NMDS.png",plot = ff, device = "png", path = "figures/", width = 12, height=8, units = "in", dpi = 400, bg = "white")


ww<-specaccum(fishdat.t,method="exact", permutations = 1000)

tidy_specaccum <- function(x) {
  data.frame(
    site = x$sites,
    richness = x$richness,
    sd = x$sd)
}
www <- tidy_specaccum(ww)

p6<-ggplot() +
  geom_line(data=yyy, aes(x=site, y=richness), linewidth=2, color="firebrick") +
  geom_linerange(data=yyy,aes(x = site, ymin = richness - 2*sd, ymax = richness + 2*sd)) +
  geom_line(data=www, aes(x=site, y=richness), linewidth=2, color="dodgerblue") +
  geom_linerange(data=www,aes(x = site, ymin = richness - 2*sd, ymax = richness + 2*sd)) +
  ylim(0, NA)+
  ylab(label = "Species Richness")+
  xlab(label="Site")+
  theme_bw()+
  theme(text = element_text(size=20))
p6
ggsave(filename = "2022SAB_COIand16S_Specaccum.png", plot = p6, device = "png", path = "figures/", width = 10, height=8, units="in",dpi = 400, bg="white")



# 2023 ESI Perley Data ----------------------------------------------------

#load sample metadata 
per23.metadata <- read.table("data/2023Perley/2023Perley-sample-metadata_ESI.tsv", sep="\t",header = T) %>% glimpse()

per23.metadata$sample.id<-gsub("-",".",per23.metadata$sample.id)

## 12S 
head(esi23.12s.perl.merge)
#Remove columns we don't need
esi23.12s.per.merge3 <-esi23.12s.per.merge2 %>% select(-c("ASV","V2","V3","V4","V5","V7","V8"))
#Group by to make some stats easier, but we should run the NMDS on the raw ASVs not grouped as species
esi23.per.smat <- esi23.12s.per.merge3 %>% group_by(V6) %>% summarise(across(everything(), sum)) %>% data.frame()



#Make a barplot of taxa
tt<-pivot_longer(esi23.per.smat, cols=starts_with("Sample"))
tt$V6<-gsub(pattern = "Clupea pallasii", replacement = "Clupea harengus", tt$V6)
tt$V6<-gsub(pattern = "Alosa fallax", replacement = "Alosa sp.", tt$V6)
tt$V6<-gsub(pattern = "Alosa pseudoharengus", replacement = "A. pseudoharengus", tt$V6)


h<- ggplot()+geom_bar(data=tt%>%filter(value>200), aes(x=V6, y=log(value)),stat="identity")+
  xlab(label = "")+
  ylab(label="12S Log(Read Count)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=55, hjust=1), text=element_text(size=14))+
  ggtitle('2023 Offshore');h

ggsave(filename = "ESI23_Perley_12S_barplot.png",plot = h, device = "png", path = "figures/2023_Perley/", width = 12, height=8, units = "in", dpi = 400, bg = "white")    

# Do NMDS of data and group by site or season

#Set up matrix for NMDS with the non-grouped data
esi23tt <- t(esi23.12s.per.merge3[,1:length(colnames(esi23.per.smat))-1])
colnames(esi23tt)<-esi23.12s.per.merge3$V6
esi23ttt<-esi23tt[rowSums(esi23tt[])>0,]

#Remove blanks - samples 22, 44, 45, and 52-56 for this dataset
esi23.mm2 <-esi23ttt[-c(15,39:40,48:49),]
nmds.esi23.fish <- metaMDS(esi23.mm2, distance = "jaccard", k=12, trymax = 200, maxit=500)

#extract nmds scores for ggplot
data.scores = as.data.frame(scores(nmds.esi23.fish)$sites)

data.scores$Sample <- rownames(data.scores)

#merge with metadata by sample id column
data.scores.merge <- inner_join(data.scores,per23.metadata, by=c("Sample"="sample.id"))


species.scores <- as.data.frame(scores(nmds.esi23.fish, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.data.per23 <- data.scores.merge %>%
  as.data.frame() %>%
  group_by(season) %>%
  slice(chull(x=NMDS1,y=NMDS2)) #for this dataset, there's no differentiation between bottom and surface
colorScales <- c("#c43b3b", "#80c43b", "#3bc4c4", "#7f3bc4")

v = ggplot() + 
  geom_polygon(data=hull.data.per23,aes(x=NMDS1,y=NMDS2,fill=season),color="black",alpha=0.30) +
  scale_fill_manual(values=c("firebrick","cornflowerblue","forestgreen"))+
  ggnewscale::new_scale_fill()+
  geom_point(data=data.scores.merge, aes(x = NMDS1, y = NMDS2, fill=station),size = 4, shape=21, colour="black")+ 
  scale_fill_manual(values=colorScales)+
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
        #legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "none", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        #legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = "NMDS1", y = "NMDS2")  + 
  geom_text(aes(x=Inf, y=Inf, vjust=48,hjust=1.2,label=paste("Stress =",round(nmds.esi23.fish$stress,3),"k =",nmds.esi23.fish$ndim)))+
  nmdstheme;v


ggsave(filename = "ESI23_Perley_12S_NMDS1_2_bySeason.png",plot = v, device = "png", path = "figures/2023_Perley/", width = 12, height=8, units = "in", dpi = 400, bg = "white")  

# Accumulation curves

esi23.12s.spec.per <-specaccum(esi23ttt,method="exact", permutations = 10000)

tidy_specaccum <- function(x) {
  data.frame(
    site = x$sites,
    richness = x$richness,
    sd = x$sd)
}
esi23.12s.spec.tidy <- tidy_specaccum(esi23.12s.spec.per)

p21<-ggplot() +
  geom_line(data=esi23.12s.spec.tidy, aes(x=site, y=richness), linewidth=2, color="firebrick") +
  geom_linerange(data=esi23.12s.spec.tidy,aes(x = site, ymin = richness - 2*sd, ymax = richness + 2*sd)) +
  ylim(0, NA)+
  ylab(label = "Species Richness")+
  xlab(label="Site")+
  theme_bw()+
  theme(text = element_text(size=20));p21

ggsave(filename = "ESI_Perley_2023_12S_SpecAccum.png",plot = p21, device = "png", path = "figures/2024CSAS/", width = 10, height=8, units = "in", dpi = 400, bg = "white")

#Alluvial plot 
#merge long data frame with metadata
jj <- left_join(tt %>% filter(!name %in% c("Sample.22","Sample.44","Sample.45","Sample.52","Sample.53")), per23.metadata,
                by=c("name"="sample.id"), 
                multiple = "all")
colnames(jj) <- c("Species", "Sample","value","watersample","conc","Depth", "Surface","Station","Season","preservative")
jj$Surface <- gsub("Deep","Bottom",jj$Surface)

#jj$Species <- str_wrap(string = jj$Species,width = 8)

filteredjj <- jj %>% filter(value>8000)

w <- ggplot(data=filteredjj, aes(axis1=Species, axis2=Season, axis3=Surface, y=log(value)))+
  scale_x_discrete(limits=c("Species","Season", "Depth"))+
  geom_alluvium(aes(fill=Station))+
  geom_stratum(alpha=0.5,width = 1/3)+  
  #geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=3)+
  ylab(label = "Log(Read count)")+
  theme_minimal()+
  theme(text=element_text(size=20));w


ggsave("Perley2023_12S_AlluvialPlot_NOLABELS.png", plot = w, device = "png", path = "figures/2023_Perley/", width = 14, height=12, dpi = 300, bg="white")



## COI 
head(esi23.coi.perl.filt)

#Remove columns we don't need
esi23.coi.per.merge3 <-esi23.coi.perl.filt %>% select(-c("ASV","Phylum","V14","Class","V16","V17","V18","V19","V20","V21","V23","V24","V26","V29"))
#Group by to make some stats easier, but we should run the NMDS on the raw ASVs not grouped as species
esi23.coi.smat <- esi23.coi.per.merge3 %>% group_by(Species) %>% summarise(across(everything(), sum)) %>% data.frame()


#Make a barplot of taxa
tt<-pivot_longer(esi23.coi.smat, cols=starts_with("Sample"))
tt$Species <- gsub("Eunoe_sp._BOLD:AAG5099","Eunoe_sp.",tt$Species)
tt$Species <- gsub("Ectyonopsis_pluridentata","Myxillidae sp.",tt$Species)
tt$Species <- gsub("Eubranchus_olivaceus","Eubranchus_sp.",tt$Species)


hh<- ggplot()+geom_bar(data=tt%>%filter(value>200,!Species=="Homo_sapiens"), aes(x=Species, y=log(value)),stat="identity")+
  xlab(label = "")+
  ylab(label="COI Log(Read Count)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=55, hjust=1), text=element_text(size=14));hh

h/hh

ggsave(filename = "ESI23_Perley_12S_COIcombined_barplot.png",plot = last_plot(), device = "png", path = "figures/2023_Perley/", width = 12, height=8, units = "in", dpi = 400, bg = "white")    

# Do NMDS of data and group by site or season

#Set up matrix for NMDS with the non-grouped data
esi23tt <- t(esi23.coi.per.merge3[,1:length(colnames(esi23.coi.per.merge3))-1])
colnames(esi23tt)<-esi23.coi.per.merge3$Species
esi23ttt<-esi23tt[rowSums(esi23tt[])>0,]

#Remove blanks - samples 22, 44, 45, and 52-56 for this dataset
esi23.mm2 <-esi23ttt[-c(15,39:40,48:52),]
nmds.esi23.inverts <- metaMDS(esi23.mm2, distance = "jaccard", k=10, trymax = 200, maxit=500)

#extract nmds scores for ggplot
data.scores = as.data.frame(scores(nmds.esi23.inverts)$sites)

data.scores$Sample <- rownames(data.scores)

#merge with metadata by sample id column
data.scores.merge <- inner_join(data.scores,per23.metadata, by=c("Sample"="sample.id"))


species.scores <- as.data.frame(scores(nmds.esi23.inverts, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.data.per23 <- data.scores.merge %>%
  as.data.frame() %>%
  group_by(season) %>%
  slice(chull(x=NMDS1,y=NMDS2)) #for this dataset, there's no differentiation between bottom and surface

vv = ggplot() + 
  geom_polygon(data=hull.data.per23,aes(x=NMDS1,y=NMDS2,fill=season),color="black",alpha=0.30) +
  scale_fill_manual(values=c("firebrick","cornflowerblue","forestgreen"))+
  ggnewscale::new_scale_fill()+
  geom_point(data=data.scores.merge, aes(x = NMDS1, y = NMDS2, fill=station),size = 4, shape=21, colour="black")+ 
  scale_fill_manual(values=colorScales)+
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2)) + 
  labs(x = "NMDS1", y = "NMDS2")  + 
  geom_text(aes(x=Inf, y=Inf, vjust=48,hjust=1.2,label=paste("Stress =",round(nmds.esi23.inverts$stress,3),"k =",nmds.esi23.inverts$ndim)))+
  nmdstheme;vv

v <- v+plot_annotation(title = "12S Fish")
vv <- vv + plot_annotation(title = "COI Eukaryotes")
v+vv+plot_annotation(title=NULL)

ggsave(filename = "ESI23_Perley_COI_12Scombined_NMDS1_2_bySeason.png",plot = last_plot(), device = "png", path = "figures/2023_Perley/", width = 12, height=8, units = "in", dpi = 400, bg = "white")  

# Accumulation curves

esi23.coi.spec.per <-specaccum(esi23ttt,method="exact", permutations = 10000)

tidy_specaccum <- function(x) {
  data.frame(
    site = x$sites,
    richness = x$richness,
    sd = x$sd)
}
esi23.coi.spec.tidy <- tidy_specaccum(esi23.coi.spec.per)

p22<-ggplot() +
  geom_line(data=esi23.coi.spec.tidy, aes(x=site, y=richness), linewidth=2, color="firebrick") +
  geom_linerange(data=esi23.coi.spec.tidy,aes(x = site, ymin = richness - 2*sd, ymax = richness + 2*sd)) +
  ylim(0, NA)+
  ylab(label = "Species Richness")+
  xlab(label="Site")+
  theme_bw()+
  theme(text = element_text(size=20));p22

ggsave(filename = "ESI_Perley_2023_COI_SpecAccum.png",plot = p22, device = "png", path = "figures/2024CSAS/", width = 10, height=8, units = "in", dpi = 400, bg = "white")

#Alluvial plot 
#merge long data frame with metadata
jj <- left_join(tt %>% filter(!name %in% c("Sample.22","Sample.44","Sample.45","Sample.52","Sample.53","Sample.54","Sample.55","Sample.56")), per23.metadata,
                by=c("name"="sample.id"), 
                multiple = "all")
colnames(jj) <- c("Species", "Sample","value","watersample","conc","Depth", "Surface","Station","Season","preservative")
jj$Surface <- gsub("Deep","Bottom",jj$Surface)

#jj$Species <- str_wrap(string = jj$Species,width = 8)

filteredjj <- jj %>% filter(value>4000)

ww <- ggplot(data=filteredjj, aes(axis1=Species, axis2=Season, axis3=Surface, y=log(value)))+
  scale_x_discrete(limits=c("Species","Season", "Depth"))+
  geom_alluvium(aes(fill=Station))+
  geom_stratum(alpha=0.5,width = 1/3)+  
  #geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=3)+
  ylab(label = "Log(Read count)")+
  theme_minimal()+
  theme(text=element_text(size=20));ww


ggsave("Perley2023_COI_AlluvialPlot_NOLABELS.png", plot = ww, device = "png", path = "figures/2023_Perley/", width = 14, height=12, dpi = 300, bg="white")

# 2023 ESI Seining Data ---------------------------------------------------

#12S 
#Remove columns we don't need from the mmm thing made in the asv_filtering.R script
mmmm<-mmm[,-c(1,3,4,5)]
#Group by to make some stats easier, but we should run the NMDS on the raw ASVs not grouped as species
esi23smat <- mmmm %>% group_by(Species) %>% summarise(across(everything(), sum)) %>% data.frame()
esi23tt <- t(esi23smat[,2:length(colnames(esi23smat))])
colnames(esi23tt)<-esi23smat$Species
esi23ttt<-esi23tt[rowSums(esi23tt[])>0,]

esi23.metadata <- read.table("data/2023Seining/seining2023-sample-metadata.tsv", sep="\t",header = T) %>% glimpse()
#don't include field blanks
groupz<-c(rep("Spring",14),rep("Summer",14), rep("Fall",15),   rep("Summer",2))
sample.sites<-c("LH","MOS","CON","MOS","GOLD","CON","CON","TAY","GOLD","LH","TAY","GOLD","LH","TAY",rep("TAY",3),rep("CON",3),rep("LH",3),rep("GOLD",3),"MOS","MOS",rep("TAY",3), rep("MOS",3),rep("GOLD",3),rep("LH",3),rep("CON",3),rep("WOLF",2))

#Make a barplot of taxa
tt<-pivot_longer(mmmm, cols=starts_with("X"))
tt$Species<-gsub(pattern = "Clupea pallasii", replacement = "Clupea harengus", tt$Species)

  p100 <- ggplot()+geom_bar(data=tt%>%filter(value>200, !Species=="Oncorhynchus keta"), aes(x=Species, y=log(value)),stat="identity")+
    xlab(label = "")+
    ylab(label="12S Log(Read Count)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust=1), text=element_text(size=14))+
    ggtitle('2023 Seining');p100

  ggsave(filename = "ESI23_Seining_12S_barplot.png",plot = last_plot(), 
         device = "png", path = "figures/2023Seining/", width = 10, height=8, 
         units = "in", dpi = 400, bg = "white")    

#Run the NMDS        
#transpose the non-grouped ASV table first
esi23.m <-t(mmmm[,2:length(colnames(mmmm))])
colnames(esi23.m) <-mmmm$Species
esi23.m<-esi23.m[rowSums(esi23.m[])>0,]
#Remove field blanks
esi23.mm <-esi23.m[-c(6,11:12,31,39,43,50,53),]
nmds.esi23.fish <- metaMDS(esi23.mm, distance = "jaccard", k=8, trymax = 200, maxit=500)
#plot(nmds.esi23.fish)


#extract nmds scores for ggplot
data.scores = as.data.frame(scores(nmds.esi23.fish)$sites)

data.scores$Sample <- rownames(data.scores)
data.scores$Season <- groupz
data.scores$Location<-sample.sites



species.scores <- as.data.frame(scores(nmds.esi23.fish, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.data <- data.scores %>%
  as.data.frame() %>%
  group_by(Season) %>%
  slice(chull(x=NMDS1,y=NMDS2))


r = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Location),color="black",alpha=0.30) + # add the hulls
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
  geom_point(size = 4, shape=21, aes(fill=Location),colour="black")+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2")  + 
  geom_text(aes(x=Inf, y=Inf, vjust=26,hjust=1.2,label=paste("Stress =",round(nmds.esi23.fish$stress,3),"k =",nmds.esi23.fish$ndim)))+
  nmdstheme;r

u = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Season),color="black",alpha=0.30) + # add the hulls
  scale_fill_manual(values=c("firebrick","cornflowerblue","forestgreen"))+
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
  geom_point(size = 4, shape=21, aes(fill=Season),colour="black")+ 
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2")  + 
  #geom_text(aes(x=Inf, y=Inf, vjust=26,hjust=1.2,label=paste("Stress =",round(nmds.esi23.fish$stress,3),"k =",nmds.esi23.fish$ndim)))+
  nmdstheme;u

u/r

ggsave(filename = "ESISeining_2023_12S_NMDS1_2_Jaccard_Site_and_Location.png",plot = last_plot(), device = "png", path = "figures/2024CSAS/", width = 10, height=8, units = "in", dpi = 400, bg = "white")

#Alluvial plot 
#merge long data frame with metadata
jj <- mmmm %>% 
      pivot_longer(cols=starts_with("X")) %>% 
      filter(!name %in% c("X2023ESI_13","X2023ESI_05","X2023ESI_19","X2023ESI_07","X499655","X499618","X499622","X499626","X499630","X499634","X499638","X499642","X499646","X499650","X499654","X499655"))

jj$name <- gsub("X","",jj$name)
jjj <- left_join(jj, esi23.metadata,
                by=c("name"="Sample.id"), 
                multiple = "all")

jjj$Species <- str_wrap(string = jjj$Species,width = 10)
jjj$site <- gsub("MH","MOS",jjj$site)
jjj$site <- gsub("OWL", "LH", jjj$site)
jjj$site <- gsub("TH", "TAY", jjj$site)


  p50 <- ggplot(data=jjj %>% filter(value>100000, !Species=="Clupea\npallasii"), aes(axis1=Species, axis2=site, axis3=season, y=log(value)))+
    scale_x_discrete(limits=c("Species","Site", "Season"))+
    geom_alluvium(aes(fill=site))+
    geom_stratum(alpha=0.5,width = 1/3)+  
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size=3)+
    ylab(label = "Log(Read count)")+
    theme_minimal()+
    theme(text=element_text(size=20));p50

ggsave("2023_Seining_12S_Alluvialplot.png",plot=p50, path="figures/2023Seining/", device = "png",
       bg="white",width = 12, height=8, dpi = 300)

#Accumulation curves

esi23.12s.spec <-specaccum(esi23.mm,method="exact", permutations = 10000)

tidy_specaccum <- function(x) {
  data.frame(
    site = x$sites,
    richness = x$richness,
    sd = x$sd)
}
esi23.12s.spec.tidy <- tidy_specaccum(esi23.12s.spec)

  p20<-ggplot() +
    geom_line(data=esi23.12s.spec.tidy, aes(x=site, y=richness), linewidth=2, color="firebrick") +
    geom_linerange(data=esi23.12s.spec.tidy,aes(x = site, ymin = richness - 2*sd, ymax = richness + 2*sd)) +
    ylim(0, NA)+
    ylab(label = "Species Richness")+
    xlab(label="Site")+
    theme_bw()+
    theme(text = element_text(size=20));p20

ggsave(filename = "ESISeining_2023_12S_SpecAccum.png",plot = p20, device = "png", path = "figures/2024CSAS/", width = 10, height=8, units = "in", dpi = 400, bg = "white")

### Make a Venn diagram of fish caught vs fish detected in eDNA
fish.catch <-read.csv("~/GitHub/easternshoreislands_aoi/data/Seining/Seining_FishMeasurements_2023.csv") %>% filter(FunctionalGroup=="BonyFish") %>% select(SpeciesName)  %>% unique()

fish.catch$SpeciesName<- gsub("Gadidae","Microgadus tomcod", fish.catch$SpeciesName)
fish.catch$SpeciesName<- gsub("Cottoidei","Myoxocephalus spp.", fish.catch$SpeciesName)
fish.catch$Spe9055ciesName<- gsub("Pholis sp.","Pholis gunnellus", fish.catch$SpeciesName)

mifish.fish$. <-gsub("Myoxocephalus scorpius","Myoxocephalus spp.", mifish.fish$.)
mifish.fish<-unique(tt$Species) %>% as.data.frame()

fish.match<-list(Seining=unique(fish.catch$SpeciesName), eDNA=mifish.fish$.)

p21<- plot(euler(fish.match), fills=list(fill=c("Seining"="#56B4E9",
                                                "eDNA"="firebrick"), alpha=0.9),
           labels=list(col="black",font=2),
           quantities=list(col="black", font=2))

p20 + p21

ggsave(filename = "ESI2023_Seining_12S_SpecAccum_andVenn.png", plot = last_plot(),device = "png", width = 12, height=8, dpi = 320, path = "figures/2023Seining/" )

#Look at alpha diversity
shan <- data.frame(diversity(esi23.mm, index="shannon"))
shan <-shan %>% rename(Shannon=diversity.esi23.mm..index....shannon..)
#we can aggregate these so we get a site level index
shan$site <-sample.sites
shan$season<-groupz
colnames(shan) <- c("Shannon","site","season")


levels=c("Spring","Summer","Fall")
ggplot()+
  geom_boxplot(data=shan,aes(x=site, y=Shannon,fill=site),alpha=0.7)+
  #geom_point(data=shan, aes(x=site,y=Shannon))+
  facet_wrap(.~factor(season, level=levels))+
  labs(x="Site",y="Shannon Diversity")+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"),
        text=element_text(size=16))


ggsave("Seining2023_ESI_ShannonDiversity_Seasonal.png",plot=last_plot(), device = "png", width=10, height=8, path = "figures/2023Seining/", dpi = 300, bg = "white")

#get species richness as n species too


df_rich_rep <- esi23.mm %>%
  rownames_to_column(var="site") %>%
  mutate(Code=substr(site, start=1, stop=3)) %>%
  group_by(Code) %>%
  summarize(avg_richness = mean(rowSums(across(where(is.numeric))>0))) %>%
  ungroup() %>%
  data.frame()

### 2023 Seining COI data - read in from 00_asv_filtering.R

head(esi23.coi.filt2)
head(esi23.metadata) #looking at the metadata again to remind myself of the object name

#remove field blanks (ESI_13, ESI_05, ESI_19, ESI_07, 499655, 499618, 499622, 499630, 499634, 499638, 499646, 499650, 499654)
esi23.coi.filt3 <- esi23.coi.filt2 %>% select(-c("X2023ESI_13","X2023ESI_19","X2023ESI_07","X2023ESI_05","X499655","X499618","X499622","X499626","X499630","X499634","X499638","X499642","X499646","X499650","X499654","X499655"))

##Create a stacked barplot per site of animal class 
#Make a barplot of taxa
cc<-pivot_longer(esi23.coi.filt3, cols=starts_with("X"))
cc$Species <- gsub("_CMC01","",cc$Species)
cc$Species <- gsub("_CMC02","",cc$Species)


  p30 <- ggplot()+geom_bar(data=cc%>%filter(value>10000), aes(x=Species, y=log(value)),stat="identity")+
    xlab(label = "")+
    ylab(label="COI Log(Read Count)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=60, hjust=1), text=element_text(size=14))+
    ggtitle('2023 Seining - COI');p30
  
  p31 <- ggplot()+geom_bar(data=cc%>%filter(!Species %in% c("Micromonas_pusilla","Acartia_hudsonica")), aes(x=Species, y=log(value)),stat="identity")+
    xlab(label = "")+
    ylab(label="COI Log(Read Count)")+
    theme_bw()+
    theme(axis.text.x = element_text(angle=60, hjust=1), text=element_text(size=14));p31
  
  p30/p31

ggsave(filename = "ESI23_Seining_COI_barplot_byAbundance.png",plot = last_plot(), device = "png", path = "figures/2023Seining/", width = 10, height=8, units = "in", dpi = 400, bg = "white") 

   # Merge long pivot table with metadata to make a stacked barplot
cc$name <- gsub("X","",cc$name)
ccc <- left_join(cc,esi23.metadata, by=c("name"="Sample.id"))

coi_summary <- ccc %>%
  filter(!site=="blank", !Class=="Mamiellophyceae") %>%
  group_by(site,Class) %>%
  summarise(total_value=sum(value),.groups="drop") %>% 
  group_by(site) %>% 
  mutate(proportion = total_value/sum(total_value))

  coi_summary$site <- gsub("OWL","LH",coi_summary$site)
  coi_summary$site <- gsub("MH","MOS",coi_summary$site)
  coi_summary$site <- gsub("TH","TAY",coi_summary$site)


  level_order=c("CON","LH","TAY","WLF","MOS","GOLD")

  colourCount = length(unique(coi_summary$Class)) #19 classes

  p99 <- ggplot(coi_summary, aes(x = factor(site,levels = level_order), y = proportion, fill = Class)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values=unname(glasbey())) +
    labs(x = "Site",
         y = "Proportion") +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal()+
    theme(text=element_text(size=16),
          axis.text.x=element_text(face="bold"));p99

ggsave(filename = "ESI23_Seining_Site_ClassBarPlots.png", plot = p99, device = "png", path = "figures/2023Seining/", width = 10, height = 8, units = "in", dpi = 300, bg = "white")

# NMDS of COI data
esi23.coi.mat <- esi23.coi.filt3 %>% 
  select(starts_with("X")) %>%
  t()

colnames(esi23.coi.mat) <-esi23.coi.filt3$Species

esi23.coi.mat<-esi23.coi.mat[rowSums(esi23.coi.mat[])>0,]
#Remove field blanks
nmds.esi23.coi <- metaMDS(esi23.coi.mat, distance = "jaccard", k=14, trymax = 200)

#extract nmds scores for ggplot
data.scores = as.data.frame(scores(nmds.esi23.coi)$sites) 


data.scores$Sample <- rownames(data.scores)
data.scores$Sample <- gsub("X","",data.scores$Sample)
data.scores2 <- left_join(data.scores,esi23.metadata, by=c("Sample"="Sample.id")
data.scores2 <- data.scores2 %>% rename(Location=site)


species.scores <- as.data.frame(scores(esi23.coi.mat, "species"))  #Using the scores function from vegan to extract the species scores and convert to a data.frame
species.scores$species <- rownames(species.scores) 

hull.data <- data.scores2 %>%
  as.data.frame() %>%
  group_by(season) %>%
  slice(chull(x=NMDS1,y=NMDS2))


  p200 = ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=season),color="black",alpha=0.30) + # add the      hulls
    scale_fill_manual(values=c("firebrick","cornflowerblue","forestgreen"),name="Season")+
    #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
    geom_point(size = 4, shape=21, aes(fill=season),colour="black")+ 
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) + 
    labs(x = "NMDS1", y = "NMDS2")  + 
    geom_text(aes(x=Inf, y=Inf, vjust=29,hjust=1.2,
                  label=paste("Stress =",
                              round(nmds.esi23.coi$stress,3),
                              "k =",nmds.esi23.coi$ndim)))+
    nmdstheme;p200

  
  hull.data <- data.scores2 %>%
    as.data.frame() %>%
    group_by(Location) %>%
    slice(chull(x=NMDS1,y=NMDS2))
  
  
  p201 = ggplot(data.scores2, aes(x = NMDS1, y = NMDS2)) + 
    geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Location),color="black",alpha=0.30) + # add the      hulls
    #scale_fill_manual(values=c("firebrick","cornflowerblue","forestgreen"))+
    #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
    geom_point(size = 4, shape=21, aes(fill=Location),colour="black")+ 
    theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
          axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
          legend.text = element_text(size = 12, face ="bold", colour ="black"), 
          legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
          axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
          legend.title = element_text(size = 14, colour = "black", face = "bold"), 
          panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
          legend.key=element_blank()) + 
    labs(x = "NMDS1", y = "NMDS2")  + 
    geom_text(aes(x=Inf, y=Inf, vjust=29,hjust=1.2,
                  label=paste("Stress =",
                              round(nmds.esi23.coi$stress,3),
                              "k =",nmds.esi23.coi$ndim)))+
    nmdstheme;p201

p200/p201

ggsave("figures/2023Seining/ESI2023_Seining_COI_NMDS1_2_Jaccard_2Panels.png", plot = last_plot(), device = "png",width = 12, height = 10,dpi = 300,bg = "white")
# Save RData --------------------------------------------------------------

  save.image(file = "data/eDNA_NMDS_and_Diversity.RData")