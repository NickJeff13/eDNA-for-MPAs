library(vegan)
library(dplyr)
library(janitor)
library(ggplot2)

#read in our data table for species accummulation curves and NMDS plots
taxtable <- read.table("data/2022Data/SAB/COI/COI_FilteredASVtable.txt", header = T)

#Make a barplot of taxa
tt<-pivot_longer(taxtable, cols=starts_with("X"))
ggplot()+geom_bar(data=tt, aes(x=V27, y=log(value)),stat="identity")+
  xlab(label = "Species")+
  ylab(label="Log(Read Count)")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))


meta <- read.table("data/2022Data/SAB/R_2022-sample-metadata_SABonly.tsv", header = T)
rownames(meta)<-meta$watersample
meta<-meta[1:79,]
head(meta)
#select some columns for this specific analysis 
#Sample key
#Aspy Bay = Samples 1-3
#Chebogue = Samples 5-7
#Lingan = Samples 8-10
#Mira = Samples 12-14
#Cheticamp = Samples 16-18
#North R = Samples 20-22
#EastBay = Samples 24-26

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


xx = ggplot(data.scores, aes(x = NMDS1, y = NMDS2)) + 
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

ggsave(filename = "2022SAB_COI_Specaccum.png", plot = p3, device = "png", path = "figures/", width = 10, height=8, units="in",dpi = 400, bg="white")

##############################################################################################
###Now do the same thing for 12S! 
######################################
fishdat <- read.table("data/2022Data/SAB/16S/16S_FilteredASVtable.txt", sep="\t",header = T)

head(fishdat)
#get rid of row 1 which is non-fish animals
fishdat <- fishdat[-1,]

fishdat.t<-t(fishdat[,2:length(colnames(fishdat))])
fishdat.t<-fishdat.t[rowSums(fishdat.t[])>0,]
colnames(fishdat.t)<-fishdat[,1]

sab.16s.nmds <-metaMDS(fishdat.t, distance="bray", k=3, trymax = 300, maxit=500)
plot(sab.16s.nmds)



data.scores.16s = as.data.frame(scores(sab.16s.nmds)$sites)

data.scores.16s$Sample <- rownames(data.scores.16s)
data.scores.16s$Depth <- groups$V2
data.scores.16s$Surface <- groups$V3


ff = ggplot(data.scores.16s, aes(x = NMDS1, y = NMDS2)) + 
  geom_point(size = 4, aes(shape = Surface, colour = Depth))+ 
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
