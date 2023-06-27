setwd("C:/Users/VANWYNGAARDENMA/Documents/eDNA/Musquash/Data/12S/dada2out_12S_Test4")
library("ggplot2")
library(RColorBrewer)
library(viridis)
library("tidyverse")


#load metadata and taxon list for ASVs
metadata <- read.delim("C:/Users/VANWYNGAARDENMA/Documents/eDNA/Musquash/Data/12S/Musquash-12S-metadata_dada2.tsv", sep="\t")
taxa <- read.csv("Musquash_12S_BlastClassifierTaxConsensus_ASV.csv")


#load reads per site and calculate total reads per ASV
readsPerSite <- read.csv("Musquash_12S_ASV_table_ReadsPerSite.csv") %>% 
  left_join(select(taxa, ASV, Consensus)) %>%
  relocate(Taxon, ASV, Consensus) %>% 
  mutate_at(c(4:90), as.numeric) %>% #convert the read data to numeric
  rowwise() %>% #ensure the rest of the command is done across each row
  mutate(TotalReads = sum(c_across(4:90)), #sum all the reads per ASV across all samples
         OnePercTotalReads = round(TotalReads*0.01,0)) %>% #calculate 1% of total reads per ASV
  relocate(Taxon, ASV, Consensus, TotalReads, OnePercTotalReads)  #reorder the columns
  
#remove any ASV counts per site that are less than 1% total reads for the ASV
readsPerSite_clean <- readsPerSite %>% 
  rowwise() %>% #ensure the rest of the command is done across each row
  mutate_at(vars(c(6:92)), funs(ifelse(.<OnePercTotalReads,0,.))) #in cols 6:92, if the value is < OnePercTotalReads for that row, replace with 0
  
#make a table that just lists the site, then which species it has in it
readsPerSite_clean$Taxon <- factor(readsPerSite_clean$Taxon)

readsPerSite_clean_long <- readsPerSite_clean %>%
  gather(Sample, Reads, X2A01S1:PCR_blank_002, factor_key=TRUE) %>% 
  relocate(Sample,Taxon,Reads,ASV) %>%
  filter(Reads !=0) %>%
  mutate(Sample = str_replace(Sample, "X","")) %>% 
  left_join(metadata)
  
write.csv(readsPerSite_clean_long, "Musquash_12S_ASV_table_ReadsPerSite_clean_long.csv")



#use bray-curtis matrix
BC_PCOACoord <- read.csv("C:/Users/VANWYNGAARDENMA/Documents/eDNA/Musquash/Data/12S/dada2out_12S_Test4/12S-core-metrics-results/bray_curtis_PCOA_coord.csv") %>% 
  left_join(metadata) %>% 
  relocate(Sample,type,location,site,zone) %>% 
  mutate_at(vars(site), factor) %>% 
  select(Sample,type,location,site,zone,P1,P2) %>% 
  mutate(site = str_replace(site, "PCR_blank_002","PCRBlank"))
  
  
SitePlot <- ggplot(BC_PCOACoord, #the dataframe that has your plot data in it
               aes(x=P1,  #column name for the x Axis data - PCA Axis1
                   y=P2)) +  #column name for the y axis data - PCA Axis2
  geom_hline(linewidth=1, colour = "black", yintercept = 0) + #horizontal line through the origin
  geom_vline(linewidth=1, colour = "black", xintercept = 0)  + #vertical line through the origin
  geom_point(aes(fill = site, 
                 shape = site), #add sample site points, fill in colour based on population
             colour ="black", #outline of points
             size=4, #size of points
             # shape = 21, #shape (google ggplot2 point shapes) - this one has both colour (the outline) and fill (the colour of the middle of the point)
             alpha = 0.7)+#transparency of the individual points
  scale_shape_manual(values = c(rep(24, times=2),rep(22, times = 3),rep(23, times = 3),rep(21, times = 17),rep(25,times=1))) +
  scale_fill_discrete() +
  # scale_colour_manual(values = palette) +
  xlab("\nPC1 (18.8%)")   +#label for the xaxis taken from the variance for PC1
  ylab("PC2 (14.1%)\n")   + #label for the yaxis taken from the variance for PC2
  theme_bw() + #standard black and white ggplot theme (makes the background white, etc)
  theme(legend.position = "right", #puts the legend at the right of the plot
        legend.title = element_blank(), #removes the title for the legend
        legend.text = element_text(size=28), #controls the text size for the legend
        legend.key = element_blank(), #removes boxes from around the individual legend elements
        axis.text = element_text(size = 28), #text size for the axis text 
        axis.title=element_text(size = 35), #text size for the axis title
        panel.border = element_rect(linewidth =0.5), #border of the plot panel
        plot.margin=unit(c(5,5,7,5), "mm"), #space around the plot, #top, right, bottom, left
        panel.grid.major=element_blank(), #remove the major grid lines in the plot background
        panel.grid.minor=element_blank()) + #remove the minor grid lines in the plot background
  NULL #this is just here so I can comment out the last line without having to remove a "+"

ggsave(SitePlot, #plot you want to save
       file = "Musquash_2022_12S_Bray-Curtis_PCoA.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width

DepthPlot <- ggplot(BC_PCOACoord, #the dataframe that has your plot data in it
                   aes(x=P1,  #column name for the x Axis data - PCA Axis1
                       y=P2)) +  #column name for the y axis data - PCA Axis2
  geom_hline(linewidth=1, colour = "black", yintercept = 0) + #horizontal line through the origin
  geom_vline(linewidth=1, colour = "black", xintercept = 0)  + #vertical line through the origin
  geom_point(aes(fill = site, 
                 shape = location), #add sample site points, fill in colour based on depth
             colour ="black", #outline of points
             size=4, #size of points
             # shape = 21, #shape (google ggplot2 point shapes) - this one has both colour (the outline) and fill (the colour of the middle of the point)
             alpha = 0.7)+#transparency of the individual points
  scale_shape_manual(values = c(23,21,22)) +
  scale_fill_discrete() +
  # scale_colour_manual(values = palette) +
  xlab("\nPC1 (18.8%)")   +#label for the xaxis taken from the variance for PC1
  ylab("PC2 (14.1%)\n")   + #label for the yaxis taken from the variance for PC2
  theme_bw() + #standard black and white ggplot theme (makes the background white, etc)
  theme(legend.position = "right", #puts the legend at the right of the plot
        legend.title = element_blank(), #removes the title for the legend
        legend.text = element_text(size=28), #controls the text size for the legend
        legend.key = element_blank(), #removes boxes from around the individual legend elements
        axis.text = element_text(size = 28), #text size for the axis text 
        axis.title=element_text(size = 35), #text size for the axis title
        panel.border = element_rect(linewidth =0.5), #border of the plot panel
        plot.margin=unit(c(5,5,7,5), "mm"), #space around the plot, #top, right, bottom, left
        panel.grid.major=element_blank(), #remove the major grid lines in the plot background
        panel.grid.minor=element_blank()) + #remove the minor grid lines in the plot background
  NULL #this is just here so I can comment out the last line without having to remove a "+"

ggsave(DepthPlot, #plot you want to save
       file = "Musquash_2022_12S_Bray-Curtis_PCoA_Depth.pdf", #filename to save the plot as
       height = 10, #height of the plot file
       width = 16,  #width of the plot file
       units = "in") #units of the height and width
