#Analyzing qPCR results and plotting on map
library(ggplot2)
library(dplyr)
library(tidyr)
library(sf)

#Read in qPCR results - both a summary of detections, and detailed results with Ct values
qpcr.sum <- read.csv("data/2024Perley/qPCR/APC0264 qPCR Results 06-Nov-24-Summary.csv", header = T) %>% 
  dplyr::select(SampleID, ABL_Sample_ID, Station, Date, Replicate,
                eDNA_Concentration, A_lupus_persample, A_lupus_overall, Wolffish_MeanCt,G_morhua_persample, 
                G_morhua_overall, Cod_MeanCt, Hippoglossus_persample, Hippoglossus_overall, Halibut_MeanCt) %>% 
  filter(Replicate==1) %>% 
  as.data.frame()

metadat <- read.csv(file = "data/2024Perley/SAB_Sampling_2024_Metadata.csv", header = T) %>% glimpse()
metadat$Lat<-as.numeric(metadat$Lat)
#merge 
qpcr.results <- left_join(qpcr.sum, y = metadat, by = "Station")
qpcr.long <- pivot_longer(data = qpcr.results, cols = c(A_lupus_overall,G_morhua_overall,Hippoglossus_overall), names_to = "Species")
#Create base map of SAB MPA
#Projections ------------
latlong <- "+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"
utm_mar <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"
CanProj <- "+proj=lcc +lat_1=49 +lat_2=77 +lat_0=63.390675 +lon_0=-91.86666666666666 +x_0=6200000 +y_0=3000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#load MPA polygons
sabshape <- read_sf("~/GitHub/stannsbank_mpa/data/Shapefiles/SAB_boundary_zones_2017.shp")%>%
  st_transform(latlong)%>%
  mutate(name="St Anns Bank")

#load benthoscape shapefile
sab_benthoscape <- read_sf("~/Github/stannsbank_mpa/data/Shapefiles/benthoscape.shp")%>%
  st_transform(latlong)%>%
  mutate(class = Assigned_c,
         class = gsub("A - ","",class), #clean up the classification labels
         class = gsub("Asp - ","",class),
         class = gsub("B - ","",class),
         class = gsub("C - ","",class),
         class = gsub("D - ","",class),
         class = gsub("E - ","",class),
         class = gsub("F - ","",class))
sab_benthoscape$class<-gsub("cobblles","cobbles",sab_benthoscape$class)

#load bathymetry 
ras <- raster::raster("~/Github/stannsbank_mpa/data/Bathymetry/sab_dem.tif") %>%
  raster::projectRaster(.,crs=latlong)
rasproj <- sp::proj4string(ras)
#basemap of Nova Scotia
novsco <- read_sf("~/GitHub/stannsbank_mpa/data/Shapefiles/NS_coastline_project_Erase1.shp")%>%st_transform(latlong)%>%
  mutate(name="Nova Scotia")%>%
  dplyr::select(name,geometry)

qPCR.plot <- ggplot()+
  geom_sf(data=novsco)+
  geom_sf(data=sab_benthoscape, aes(fill=class))+
  geom_sf(data=sabshape,fill=NA,linewidth=1,colour="black")+
  #coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.6,46.5), expand=F)+
  geom_point(data=filter(qpcr.long,!is.na(Long)), aes(x=Long, y=Lat,shape=value),size=3)+
  scale_shape_manual(values=c("Detected"=16, "Not Detected"=4, "Inconclusive"=10))+
  facet_grid(.~Species)+
  theme_bw()+
  #scale_fill_viridis(option="A",direction = -1,na.value="transparent")+
  labs(x="",y="",fill="Benthoscape")+
  theme(legend.position = "right",
        text=element_text(size=16), 
        legend.key.width = unit(1.5, "cm"),
        strip.background = element_rect(fill="white"));qPCR.plot

ggsave(filename = "SAB_qPCR_SummaryPlot.png", plot = qPCR.plot, device = "png", path = "figures/2024_StAnnsBank", create.dir = TRUE, width = 18, height=14, dpi = 300)


#Now a map with actual Ct values and points based on these values
ctmean.long <- pivot_longer(data = qpcr.results, cols = c(Wolffish_MeanCt,Cod_MeanCt, Halibut_MeanCt), names_to = "Species")

ctmean.long$Species <- gsub(pattern = "_MeanCt", replacement = "", ctmean.long$Species)
ctmean.long$Species <- gsub(pattern = "Cod", replacement = "Atlantic cod", ctmean.long$Species)
ctmean.long$Species <- gsub(pattern = "Wolffish", replacement = "Atlantic wolffish", ctmean.long$Species)

#plotting 1/value because a higher Ct means lower eDNA presence, and lower Ct means more eDNA
ctmean.plot <- ggplot()+
  geom_sf(data=novsco)+
  geom_sf(data=sab_benthoscape, aes(fill=class))+
  geom_sf(data=sabshape,fill=NA,linewidth=1,colour="black")+
  #coord_sf(xlim=plotlims[c(1,3)],ylim=plotlims[c(2,4)])+
  coord_sf(xlim=c(-60, -58.2), ylim=c(45.6,46.5), expand=F)+
  geom_point(data=filter(ctmean.long,!is.na(Long)), aes(x=Long, y=Lat, size=1/value))+
  #scale_shape_manual(values=c("Detected"=16, "Not Detected"=4, "Inconclusive"=10))+
  scale_size(guide="none")+
  facet_grid(.~Species)+
  theme_bw()+
  #scale_fill_viridis(option="A",direction = -1,na.value="transparent")+
  labs(x="",y="",fill="Benthoscape")+
  theme(legend.position = "bottom",
        text=element_text(size=14), 
        legend.key.width = unit(1.5, "cm"),
        strip.background = element_rect(fill="white"));ctmean.plot


ggsave(filename = "SAB_2024_qPCR_MeanCtMap.png", plot = ctmean.plot, device = "png", path = "figures/2024_StAnnsBank/", width = 12, height=8, dpi=300)
