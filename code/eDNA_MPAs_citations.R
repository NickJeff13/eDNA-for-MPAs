#Code to create a plot showing trends in publication

#Data from Web of Science Search on March 26 2025

#Search Terms
  #  (("eDNA" OR "environmental DNA") AND ("MPA" OR "Marine Protected Area" OR "Marine Protected Areas"))
  # (("eDNA" OR "environmental DNA") AND ("metabarcoding"))
  # (("eDNA" OR "environmental DNA") AND ("metabarcoding") AND ("Marine"))

#load libraries
library(tidyverse)

#file paths
wos_paths <- dir("data/CitationSearch/",full.names = TRUE)

#What (Julien) day of the year is it
doy <- 31+28+25

wos_results <- NULL
for(i in 1:3){
  
  search_term <- gsub("data/CitationSearch/","",wos_paths[i])%>%
                 gsub(".txt","",.)%>%
                 gsub("_"," ",.)
  
  wos_results <- read.table(wos_paths[i],header=TRUE)%>%
    rename(year=1,count=2,percent=3)%>%
    mutate(search=search_term)%>%
    dplyr::select(year,count,percent,search)%>%
    arrange(-year)%>%
    mutate(count = ifelse(year==2025,count*365/doy,count))%>%
    rbind(.,wos_results)
  
}


citation_plot <- ggplot(data=wos_results%>%filter(year>2011,year<2025),aes(x=year,y=count,group=search,fill=search))+
  geom_line()+
  geom_point(shape=21,col="black",size=3)+
  theme_bw()+
  labs(fill="WoS search")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.2,0.9),
        legend.background = element_blank())

ggsave("figures/eDNA_citation_plot_March2025.png",citation_plot,height = 6,width = 6,units="in",dpi=300)


