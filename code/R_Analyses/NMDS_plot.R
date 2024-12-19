# Code for running species accumulation curves, NMDS, and other diversity metrics on a whole bunch of eDNA metabarcoding data including species richness plots, NMDS plots, and species accummulation plots. Written by Nick.Jeffery@dfo-mpo.gc.ca 

# Load libraries ----------------------------------------------------------

library(vegan)
library(dplyr)
library(janitor)
library(ggplot2)
library(tidyverse)
library(ggalluvial) #make alluvial plots
library(patchwork) #stick plots together
library(eulerr) #for Venn diagrams


# NMDS theme for all plots ------------------------------------------------

nmdstheme <- theme_bw()+
  theme(#axis.text.x=element_blank(),
    #axis.ticks.x=element_blank(),
    #axis.text.y=element_blank(),
    #axis.ticks.y=element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text=element_text(size=18))

# Start loading the data

# 2021 ESI Perley Data ----------------------------------------------------

#### metadata
esi21meta <-read.table("data/2021Data/metadata/2021-sample-metadata_ESIonly.tsv",header = T, sep = "\t")
esi21meta$sample.id<-gsub("-",".",esi21meta$sample.id)

#### 12S
glimpse(esi12s.filt)
#esi12smat <- esi12s.filt %>% group_by(Species) %>% summarise(across(everything(), sum)) %>% data.frame()
esi12tt <- t(esi12s.filt.fish[,2:length(colnames(esi12s.filt.fish))])
colnames(esi12tt)<-esi12s.filt.fish$Species
esi12ttt<-esi12tt[rowSums(esi12tt[])>0,]

##### Make a barplot of taxa
tt<-pivot_longer(esi12s.filt.fish, cols=starts_with("Sample"))
a<- ggplot()+geom_bar(data=tt%>%filter(value>2000), aes(x=Species, y=log(value)),stat="identity")+
  xlab(label = "")+
  ylab(label="12S Log(Read Count)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1), text=element_text(size=14))+
  ggtitle('2021');a

#16S
glimpse(esi16s.filt)
spec.mat<-as.data.frame(t(esi16s.filt[,2:101]))
colnames(spec.mat)<-esi16s.filt$species
spec.mat<-spec.mat[rowSums(spec.mat)>0,]

#Make a barplot of taxa
tt<-pivot_longer(esi16s.filt, cols=starts_with("Sample"))
b <- ggplot()+geom_bar(data=tt%>%filter(value>3000), aes(x=species, y=log(value)),stat="identity")+
  xlab(label = "Species")+
  ylab(label="16S Log(Read Count)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1), text=element_text(size=14));b

#plot with patchwork
a/b

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

#plot it up
ggplot(data.scores.metadat %>% filter(surface !="BLANK"), aes(x = NMDS1, y = NMDS2)) + 
  #geom_text(data=species.scores,aes(x=NMDS1,y=NMDS2,label=species), alpha=0.5)+
  geom_jitter(size = 4, aes(fill = depth, shape = surface), 
              width = 10, height = 10) +
  scale_fill_viridis_c(option="viridis") + # Continuous color scale for depth
  scale_shape_manual(values = c(21, 24)) +         # Discrete shapes (e.g., circle and triangle)
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 14), 
        legend.text = element_text(size = 12, face ="bold", colour ="black"), 
        legend.position = "right", axis.title.y = element_text(face = "bold", size = 14), 
        axis.title.x = element_text(face = "bold", size = 14, colour = "black"), 
        legend.title = element_text(size = 14, colour = "black", face = "bold"), 
        panel.background = element_blank(), panel.border = element_rect(colour = "black", fill = NA, size = 1.2),
        legend.key=element_blank()) + 
  labs(x = "NMDS1", y = "NMDS2")  + 
  geom_text(aes(x=Inf, y=Inf, vjust=41,hjust=1.1,label=paste("Stress =",round(esi12s.nmds.jac$stress,3),"k =",esi12s.nmds.jac$ndim)))

ggsave(filename = "ESI21_16S_NMDS1_NMDS1_Jaccard_FishOnly.png",plot = last_plot(), device = "png", path = "figures/2021Results/", width = 10, height=8, dpi = 320)


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



spec.mat<-as.data.frame(t(esi22.12s.merge %>% select(contains("Sample.")) %>% as.data.frame()))
colnames(spec.mat)<-esi22.12s.merge$V6
spec.mat<-spec.mat[rowSums(spec.mat)>0,]

#Make a barplot of taxa
tt<-pivot_longer(esi22.12s.merge, cols=starts_with("Sample"))


f <- ggplot()+geom_bar(data=tt%>%filter(value>2000), aes(x=V6, y=log(value)),stat="identity")+
  xlab(label = "Species")+
  ylab(label="12S Log(Read Count)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1), text=element_text(size=14))+
  ggtitle('2022');f


###### 16S Data ##
head(esi22.16s.merge)
esi22.16s.merge$V6 <- gsub("Sebastes mentella", "Sebastes sp.", esi22.16s.merge$V6)
esi22.16s.merge$V6 <- gsub("Gadus macrocephalus", "Gadus morhua", esi22.16s.merge$V6)
esi22.16s.merge$V6 <- gsub("Pholis laeta", "Pholis gunnellus", esi22.16s.merge$V6)
esi22.16s.merge$V6 <- gsub("Platichthys environmental sample", "Platichthys flesus", esi22.16s.merge$V6)
esi22.16s.merge$V6 <- gsub("Pollachius virens", "Pollachius pollachius", esi22.16s.merge$V6)
esi22.16s.merge$V6 <- gsub("Pleuronectes platessa", "Pleuronectinae", esi22.16s.merge$V6)


tt<-pivot_longer(esi22.16s.merge, cols=starts_with("Sample"))

g <- ggplot()+geom_bar(data=tt%>%filter(value>2000, !V6 %in% c("Platichthys flesus","Myzopsetta punctatissima")), aes(x=V6, y=log(value)),stat="identity")+
  xlab(label = "Species")+
  ylab(label="16S Log(Read Count)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1), text=element_text(size=14));g

#plot with patchwork
f/g

ggsave(filename = "ESI2022_FishBarplots_Combined.png",plot = last_plot(), device = "png", path = "figures/2022Results//", width = 10, height=8, dpi = 300, bg = "white")

#Run the NMDS on the 16S data


###### LerayXT COI Data ##
head(esi22.coi.merge)



# 2022 SAB Perley Data ----------------------------------------------------


#read in our data table for species accummulation curves and NMDS plots
taxtable <- read.table("data/2022Data/SAB/COI/COI_FilteredASVtable.txt", header = T)

#Make a barplot of taxa
tt<-pivot_longer(taxtable, cols=starts_with("X"))
ggplot()+geom_bar(data=tt, aes(x=V27, y=log(value)),stat="identity")+
  xlab(label = "Species")+
  ylab(label="Log(Read Count)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1, vjust=0.2), text=element_text(size=14))

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


xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
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
  geom_text(aes(x=Inf, y=Inf, vjust=65,hjust=1.2,label=paste("Stress =",round(sab.coi.nmds$stress,3),"k =",sab.coi.nmds$ndim)))+  scale_colour_continuous(trans="reverse")

xx


ggsave(filename = "SAB2022_COI_NMDS.png",plot = xx, device = "png", path = "figures/", width = 10, height=8, units = "in", dpi = 400, bg = "white")

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

p3<-ggplot() +
  geom_line(data=yyy, aes(x=site, y=richness), linewidth=2, color="firebrick") +
  geom_linerange(data=yyy,aes(x = site, ymin = richness - 2*sd, ymax = richness + 2*sd)) +
  ylim(0, NA)+
  ylab(label = "Species Richness")+
  xlab(label="Site")+
  theme_bw()+
  theme(text = element_text(size=20))
p3
ggsave(filename = "2022SAB_COI_Specaccum.png", plot = p3, device = "png", path = "figures/", width = 10, height=8, units="in",dpi = 400, bg="white")

###Now do the same thing for 12S! 
######################################
fishdat <- read.table("data/2022Data/SAB/16S/16S_FilteredASVtable.txt", sep="\t",header = T)
s.groups <- read.csv("data/2022Data/SAB/SAB16S_GroupData.csv")
head(fishdat)
#get rid of row 1 which is non-fish animals
fishdat <- fishdat[-c(1,16),]

mm<-pivot_longer(fishdat, cols=starts_with("X"))
ggplot()+geom_bar(data=mm, aes(x=Species, y=log(value)),stat="identity")+
  xlab(label = "Species")+
  ylab(label="log(Read Count)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1, vjust=1), text=element_text(size=18))

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
tt<-pivot_longer(esi23smat, cols=starts_with("X"))
tt$Species<-gsub(pattern = "Clupea pallasii", replacement = "Clupea harengus", tt$Species)

n<- ggplot()+geom_bar(data=tt%>%filter(value>200, !Species=="Oncorhynchus keta"), aes(x=Species, y=log(value)),stat="identity")+
  xlab(label = "")+
  ylab(label="12S Log(Read Count)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust=1), text=element_text(size=14))+
  ggtitle('2023 Seining');n
 
ggsave(filename = "ESI23_Seining_12S_barplot.png",plot = last_plot(), device = "png", path = "figures/2023Seining/", width = 10, height=8, units = "in", dpi = 400, bg = "white")    

#Run the NMDS        
#transpose the non-grouped ASV table first
esi23.m <-t(mmmm[,2:length(colnames(mmmm))])
colnames(esi23.m) <-mmmm$Species
esi23.m<-esi23.m[rowSums(esi23.m[])>0,]
#Remove field blanks
esi23.mm <-esi23.m[-c(6,11:12,31,39,43,50,53),]
nmds.esi23.fish <- metaMDS(esi23.mm, distance = "bray", k=10, trymax = 200, maxit=500)
plot(nmds.esi23.fish)


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
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Location, color=Location),alpha=0.30) + # add the hulls
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
  geom_text(aes(x=Inf, y=Inf, vjust=37,hjust=1.2,label=paste("Stress =",round(nmds.esi23.fish$stress,3),"k =",nmds.esi23.fish$ndim)))+
  nmdstheme;r

u = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
  geom_polygon(data=hull.data,aes(x=NMDS1,y=NMDS2,fill=Season, color=Season),alpha=0.30) + # add the hulls
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
  geom_text(aes(x=Inf, y=Inf, vjust=37,hjust=1.2,label=paste("Stress =",round(nmds.esi23.fish$stress,3),"k =",nmds.esi23.fish$ndim)))+
  nmdstheme;u

u/r

ggsave(filename = "ESISeining_2023_12S_NMDS1_2_Bray_bySite.png",plot = r, device = "png", path = "figures/2024CSAS/", width = 10, height=8, units = "in", dpi = 400, bg = "white")

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

levels=c("Spring","Summer","Fall")
ggplot()+
  geom_boxplot(data=shan,aes(x=site, y=Shannon,fill=site,alpha=0.7))+
  #geom_point(data=shan, aes(x=site,y=Shannon))+
  facet_wrap(.~factor(season, level=levels))+
  theme_bw()+
  theme(strip.background = element_rect(fill="white"))
#get species richness as n species too
rich_df<-as.data.frame(commat2)

df_rich_rep <- rich_df %>%
  rownames_to_column(var="site") %>%
  mutate(Code=substr(site, start=1, stop=3)) %>%
  group_by(Code) %>%
  summarize(avg_richness = mean(rowSums(across(where(is.numeric))>0))) %>%
  ungroup() %>%
  data.frame()

### 2023 Seining COI data - read in from 00_asv_filtering.R
head(esi23.coi.filt)

##Create a stacked barplot per site of animal class 
esi23.coi.long <- esi23.coi.filt %>% 
  gather(Site_Rep, Count,-c(OTU.ID, Phylum, Class, Species, V24,V26,V29)) %>% 
  separate(Site_Rep, into=c("Site", "Replicate"), sep="\\.")

spec_summary<- tax_long %>%
  group_by(Site, Class) %>% 
  summarise(TotalCount =sum(Count)) %>%
  ungroup()

#could probably combine this with the above piping 
spec_summary <- spec_summary %>% 
  group_by(Site) %>%
  mutate(RelativeCount = TotalCount/sum(TotalCount)) %>%
  ungroup() %>%
  mutate(Site=str_replace(Site, "WRK02","WRK"))

level_order <- c("CHB","RSE","FAI","FRK","WRK","CON","PLS","CAB","TAW","TAE","MOS","ESB","MIR","NOR","LIN","CHT","ASB")
#colour Site names by their colours in the map
axis_pal <- c(rep("#d1495b",3), rep("#edae49",8), rep("#00798c",6))
colourCount = length(unique(spec_summary$Class))
getPalette = colorRampPalette(brewer.pal(9, "Set1"))

ggplot(spec_summary, aes(x = factor(Site, level=level_order), y = RelativeCount, fill = Class)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values=getPalette(colourCount)) +
  labs(x = "Site",
       y = "Relative Count") +
  scale_y_continuous(labels = scales::percent_format()) +
  theme_minimal()+
  theme(text=element_text(size=16),
        axis.text.x=element_text(colour=axis_pal,face="bold"))

ggsave(filename = "Site_ClassBarPlots.png", plot = last_plot(),device = "png", path = "output/", width = 10, height = 8, units = "in", dpi = 300)


# 2023 ESI Perley Data ----------------------------------------------------
#12S 
head(esi23.12s.perl.merge)
#Remove columns we don't need
esi23.12s.per.merge3 <-esi23.12s.per.merge2 %>% select(-c("ASV","V2","V3","V4","V5","V7","V8"))
#Group by to make some stats easier, but we should run the NMDS on the raw ASVs not grouped as species
esi23.per.smat <- esi23.12s.per.merge3 %>% group_by(V6) %>% summarise(across(everything(), sum)) %>% data.frame()
esi23tt <- t(esi23smat[,2:length(colnames(esi23smat))])
colnames(esi23tt)<-esi23smat$Species
esi23ttt<-esi23tt[rowSums(esi23tt[])>0,]

per23.metadata <- read.table("data/2023Perley/ESI/", sep="\t",header = T) %>% glimpse()
#don't include field blanks
groupz<-c(rep("Spring",14),rep("Summer",14), rep("Fall",15),   rep("Summer",2))
sample.sites<-c("LH","MOS","CON","MOS","GOLD","CON","CON","TAY","GOLD","LH","TAY","GOLD","LH","TAY",rep("TAY",3),rep("CON",3),rep("LH",3),rep("GOLD",3),"MOS","MOS",rep("TAY",3), rep("MOS",3),rep("GOLD",3),rep("LH",3),rep("CON",3),rep("WOLF",2))

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



#Save RData
save.image(file = "data/eDNA_NMDS_and_Diversity.RData")
