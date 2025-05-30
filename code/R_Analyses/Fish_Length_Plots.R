library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

#Plotting simple histograms of fish length distributions 
fish.x <-read.csv("data/2019/2019_FishLengthMeasurements.csv", header = T)
colnames(fish.x)<-gsub("Total_Length..mm.", "length", colnames(fish.x))

fish.xx <- fish.x %>% group_by(Species) %>%
  summarise(occurrences=n()) %>% 
  filter(occurrences > 30) %>% 
  inner_join(fish.x, by="Species") %>%
  data.frame()
  
fish.xx$Species<-str_wrap(fish.xx$Species, width=8)
fish.xx$Site <- gsub("WRK02","WRK",fish.xx$Site)
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
  geom_histogram(binwidth = 8, color="black", alpha=0.7)+
  facet_grid(cols=vars(Species),scales="free_x")+
  labs(x="Fish Total Length (mm)", y="Frequency")+
  theme_minimal()+
  theme(strip.text.y=element_text(size=10, face="bold"),
        strip.text.x=element_text(size=14, face="bold", angle=45),
        axis.text = element_text(size=12),
        legend.position = "bottom", panel.spacing=unit(1, "lines"))


ggsave("2019_FishLength_Histogram_nonfacet.png",plot=last_plot(),width=16, height=10, dpi=300, path = "figures/", bg = "white")


###2023 data###

fish.y <-read.csv("~/GitHub/easternshoreislands_aoi/data/Seining/Seining_FishMeasurements_2023.csv") %>% glimpse()

fish.yy <- fish.y %>% 
  select(Site, Season, SpeciesName, FunctionalGroup, Total_Length..mm.) %>% 
  filter(FunctionalGroup=="BonyFish") %>%
  data.frame()

fish.yyy <-fish.yy %>% group_by(SpeciesName) %>%
  summarise(occurences=n()) %>%
  filter(occurences >5) %>%
  inner_join(fish.yy, by="SpeciesName") %>%
  data.frame()

#Wrap species names to fit the ggplot better
fish.yyy$SpeciesName<-str_wrap(fish.yyy$SpeciesName, width=1)
fish.yyy$SpeciesName<-str_replace_all(fish.yyy$SpeciesName, "Pseudopleuronectes", "P.")
#fish.yyy$SpeciesName<-str_replace_all(fish.yyy$SpeciesName, "Tautogolabrus", "T.")
fish.yyy$Site<-str_replace_all(fish.yyy$Site, "Conrod", "Conrods Beach")

ggplot(data=fish.yyy, aes(as.numeric(Total_Length..mm.),fill=Site))+
  geom_histogram(binwidth = 5, color="black", alpha=0.7)+
  facet_grid(rows = vars(Site), cols=vars(SpeciesName),scales=c("free"))+
 scale_fill_manual(values=c("#FC4E07","violet", "#293352", "#E69F00","#56B4E9","#009E73"))+
  labs(x="Fish Length (mm)", y="Frequency")+
  theme_minimal()+
  theme(strip.text.x=element_text(size=12, face="bold", angle=45),
        strip.text.y=element_text(size=13, face="bold", angle=0),
        axis.text = element_text(size=14),
        legend.position = "none")

ggsave("2023_FishLength_Histogram_bySite.png",plot=last_plot(), width=14, height=10, dpi=300, path="figures/2024CSAS/", bg="white")

fish.yyy$Season<-factor(fish.yyy$Season, levels=c("Spring","Summer","Fall"))

#filter to only Moosehead and TH for seasonal plot
ggplot(data=fish.yyy %>% filter(Site==c("Taylor Head","Moosehead")), aes(as.numeric(Total_Length..mm.),fill=Site))+
  geom_histogram(binwidth=5, color="black", alpha=0.8)+
  facet_grid(rows = vars(Season), cols=vars(SpeciesName),scales="free")+
  #scale_x_continuous(breaks = scales::pretty_breaks(n = 4)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  labs(x="Fish Length (mm)", y="Frequency")+
  theme_minimal()+
  theme(strip.text.x=element_text(size=11, face="bold", angle=45),
        strip.text.y=element_text(size=14, face="bold"),
        axis.text = element_text(size=14),
        text=element_text(size=20),
        legend.position = "bottom",
        panel.spacing=unit(1, "lines"))

ggsave("2023_FishLength_Histogram_FacetBySeason.png",plot=last_plot(),width=18, height=14, dpi=300, path = "figures/2024CSAS/", bg = "white")


###2024 data###

fish.u <-read.csv("~/GitHub/easternshoreislands_aoi/data/Seining/Seining_FishMeasurements_2024.csv") %>% glimpse()

fish.uu <- fish.u %>% 
  select(Site, Season, SpeciesName, FunctionalGroup, Total_Length..mm.) %>% 
  filter(FunctionalGroup=="BonyFishes", Total_Length..mm.>0) %>%
  data.frame()

fish.uuu <-fish.uu %>% group_by(SpeciesName) %>%
  summarise(occurences=n()) %>%
  filter(occurences >5) %>%
  inner_join(fish.uu, by="SpeciesName") %>%
  data.frame()

#Wrap species names to fit the ggplot better
fish.uuu$SpeciesName<-str_wrap(fish.uuu$SpeciesName, width=1)
fish.uuu$SpeciesName<-str_replace_all(fish.uuu$SpeciesName, "Pseudopleuronectes", "P.")
fish.uuu$SpeciesName<-str_replace_all(fish.uuu$SpeciesName, "Tautogolabrus", "T.")
fish.uuu$Site<-str_replace_all(fish.uuu$Site, "Owls Head", "Little Harbour")

ggplot(data=fish.uuu, aes(as.numeric(Total_Length..mm.),fill=Site))+
  geom_histogram(binwidth = 10, color="black", alpha=0.7)+
  facet_grid(rows = vars(Site), cols=vars(SpeciesName),scales="free")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_fill_manual(values=c("#FC4E07","violet", "#293352", "#E69F00","#56B4E9"))+
  labs(x="Fish Length (mm)", y="Frequency")+
  theme_minimal()+
  theme(strip.text.x=element_text(size=12, face="bold",angle=45),
        strip.text.y=element_text(size=13, face="bold",angle=0),
        axis.text = element_text(size=14),
        legend.position = "none")

ggsave("2024_FishLength_Histogram_bySite.png",plot=last_plot(), width=14, height=10, dpi=300, path="figures/2024CSAS/", bg="white")

fish.uuu$Season<-factor(fish.uuu$Season, levels=c("Spring","Summer","Fall"))

#filter to just Taylor Head and Moosehead for seasonal comparison 
ggplot(data=fish.uuu %>% filter(Site==c("Taylor Head","Moosehead")), aes(as.numeric(Total_Length..mm.),fill=Site))+
  geom_histogram(binwidth=5, color="black", alpha=0.8)+
  facet_grid(rows = vars(Season), cols=vars(SpeciesName),scales="free")+
  scale_x_continuous(breaks = scales::pretty_breaks(n = 3)) +
  scale_fill_manual(values=c("#E69F00", "#56B4E9"))+
  labs(x="Fish Length (mm)", y="Frequency")+
  theme_minimal()+
  theme(strip.text.x=element_text(size=11, face="bold", angle=45),
        strip.text.y=element_text(size=14, face="bold"),
        axis.text = element_text(size=14),
        text=element_text(size=20),
        legend.position = "bottom",
        panel.spacing=unit(1, "lines"))

ggsave("2024_FishLength_Histogram_FacetBySeason.png",plot=last_plot(),width=18, height=14, dpi=300, path = "figures/2024CSAS/", bg = "white")

