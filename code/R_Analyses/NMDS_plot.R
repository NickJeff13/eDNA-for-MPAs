library(vegan)
library(dplyr)
library(janitor)
library(ggplot2)

#read in our data table for species accummulation curves and NMDS plots
taxtable <- read.table("data/COI_TaxonomyPerSite_2023.txt", header = T)

#select some columns for this specific analysis 
#Sample key
#Aspy Bay = Samples 1-3
#Chebogue = Samples 5-7
#Lingan = Samples 8-10
#Mira = Samples 12-14
#Cheticamp = Samples 16-18
#North R = Samples 20-22
#EastBay = Samples 24-26

commat <- taxtable %>% select(V26, Sample.1,  Sample.2, Sample.3,  Sample.5, Sample.6, Sample.7, Sample.8, Sample.9,Sample.10, Sample.12, Sample.13, Sample.14, Sample.16, Sample.17, Sample.18, Sample.20, Sample.21, Sample.22, Sample.24, Sample.25, Sample.26) #removed samples 4, 11, 15, 19, 23, 27 as these are the negatives

#Transpose the table for vegan
commat2<-t(commat[,2:length(colnames(commat))])
colnames(commat2)<-commat[,1]

trask.nmds <-metaMDS(commat2, distance="bray", k=5, trymax = 200)
plot(trask.nmds) #this is not very informative without labels!

ordiplot(trask.nmds, type='n')
orditorp(trask.nmds, display="species", col="red", air=0.01)
orditorp(trask.nmds, display="sites", cex=1, air=0.01)

