library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

#Plotting simple histograms of fish length distributions 
fish.x <-read.csv("data/2019_FishLengthMeasurements.csv", header = T)
colnames(fish.x)<-gsub("Total_Length..mm.", "length", colnames(fish.x))

fish.xx <- fish.x %>% group_by(Species) %>%
  summarise(occurrences=n()) %>% 
  filter(occurrences > 40) %>% 
  inner_join(fish.x, by="Species") %>%
  data.frame()
  
fish.xx$Species<-str_wrap(fish.xx$Species, width=8)

ggplot(data=fish.xx, aes(as.numeric(length),fill=Species))+
  geom_histogram(binwidth = 5, color="black", alpha=0.7)+
  facet_grid(rows = vars(Site), cols=vars(Species),scales="free_y")+
  labs(x="Fish Length (mm)", y="Frequency")+
  theme_minimal()+
  theme(strip.text.x=element_text(size=10, face="bold"),
        strip.text.y=element_text(size=10, face="bold"),
        axis.text = element_text(size=12),
        legend.position = "none")

ggplot(data=fish.xx, aes(as.numeric(length),fill=Site))+
  geom_histogram(binwidth = 5, color="black", alpha=0.7)+
  facet_grid(cols=vars(Species),scales="free_x")+
  labs(x="Fish Length (mm)", y="Frequency")+
  theme_minimal()+
  theme(strip.text.y=element_text(size=10, face="bold"),
        axis.text = element_text(size=12),
        legend.position = "right")


ggsave("2019_FishLength_Histogram.png",plot=last_plot(),width=10, height=8, dpi=300, path = "figures/")
