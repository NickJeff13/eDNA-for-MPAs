#Make a map of eDNA sites in the ESI including our Perley work,and Kira's kelp sites
library(sf)
library(raster)
library(ggplot2)
library(dplyr)

latlong <- "+init=epsg:4326 +proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"
utm <- "+proj=utm +zone=20 +datum=NAD83 +units=km +no_defs +ellps=GRS80 +towgs84=0,0,0"

#Kira sites
sites <- data.frame(sites=c("Barren Is","Borgles Is","Goose Is", "Hallibut Is", "Harbour Is", "High Is", "Laybold Is", "Phoenix Is", "Salisbury Is", "Sober Is"), lat=c(44.69, 44.767, 44.939, 44.892, 44.87, 44.8859, 44.687, 44.782, 44.822, 44.825), long=c(-62.969, -62.723, -62.053, -62.195, -62.334, -62.3234, -62.847, -62.622, -62.516, -62.445))

#Perley transects


novsco <-read_sf("../Courtney-Trask-Honours/data/NS_coastline_project_Erase1.shp") %>% 
  st_transform(latlong) %>%
  mutate(name="Nova Scotia")%>%
  dplyr::select(name, geometry)

esi <-read_sf("../Courtney-Trask-Honours/data//EasternShoreIslands_networksite.shp")%>%
  st_transform(latlong)%>%
  mutate(name="ESI")%>%
  dplyr::select(name,geometry)


focal_bound <- esi%>%
  st_bbox()%>%
  st_as_sfc()%>%
  st_transform(utm)%>%
  st_buffer(4)%>% 
  st_transform(latlong)%>%
  st_bbox()

ggplot()+
  geom_sf(data=novsco, fill="grey")+
  geom_sf(data=esi, fill=NA, linewidth=1.25, colour="black")+
  coord_sf(expand=0, xlim=focal_bound[c(1,3)], ylim=focal_bound[c(2,4)])+
  xlab(label = "Longitude")+
  ylab(label="Latitude")+
  geom_point(data=sites, aes(x=long, y=lat), shape=21, color="black",fill="red", size=4)+
  theme_bw()+
  theme(axis.text.y = element_text(colour = "black", size = 12, face = "bold"), 
        axis.text.x = element_text(colour = "black", face = "bold", size = 12))

ggsave(filename = "Kira_eDNA_Sites.png",plot = last_plot(), width = 10, height=8, dpi=300, path = "figures/2024CSAS/")
