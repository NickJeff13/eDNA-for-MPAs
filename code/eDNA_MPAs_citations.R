#Code to create a plot showing trends in publication

#Data from Web of Science Search on March 26 2025

#Search Terms
  #  (("eDNA" OR "environmental DNA") AND ("MPA" OR "Marine Protected Area" OR "Marine Protected Areas"))
  # (("eDNA" OR "environmental DNA") AND ("metabarcoding"))
  # (("eDNA" OR "environmental DNA") AND ("metabarcoding") AND ("Marine"))

#load libraries
library(tidyverse)
library(tm)
library(wordcloud2)
library(htmlwidgets)
library(webshot)

#file paths
wos_paths <- dir("data/CitationSearch/",full.names = TRUE)
wos_paths <- wos_paths[!grepl("savedrecs",wos_paths)] #this is where the abstract information comes form 


#Publication per year -------
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
    mutate(count_org =count,
           count = ifelse(year==2025,count*365/doy,count))%>%
    rbind(.,wos_results)
  
}


citation_plot <- ggplot(data=wos_results%>%filter(year>2011,year<2025),aes(x=year,y=count,group=search,fill=search))+
  geom_line()+
  geom_point(shape=21,col="black",size=3)+
  theme_bw()+
  labs(fill="WoS search",
       x="",
       y="Citation count")+
  theme(legend.position = "inside",
        legend.position.inside = c(0.2,0.9),
        legend.background = element_blank())+
  scale_x_continuous(breaks=seq(2012,2024,2))

ggsave("figures/eDNA_citation_plot_March2025.png",citation_plot,height = 6,width = 6,units="in",dpi=300)

##word cloud ---------

data <- read.csv("data/CitationSearch/savedrecs.csv", stringsAsFactors = FALSE)

text <- paste(data$Abstract, collapse = " ")
text <- gsub("Marine Protected Areas", "Marine_Protected_Area", text, ignore.case = TRUE)
text <- gsub("Marine Protected Area", "Marine_Protected_Area", text, ignore.case = TRUE)
text <- gsub("Marine Protected", "Marine_Protected_Area", text, ignore.case = TRUE)
text <- gsub("Sharks","Shark",text,ignore.case = TRUE)
text <- gsub("Communities","Community",text,ignore.case = TRUE)

# Create corpus and clean text
corpus <- Corpus(VectorSource(text))
corpus <- tm_map(corpus, content_transformer(tolower))
corpus <- tm_map(corpus, removePunctuation)
corpus <- tm_map(corpus, removeWords, stopwords("english"))

# Create Term Document Matrix
tdm <- TermDocumentMatrix(corpus)
m <- as.matrix(tdm)
word_freqs <- sort(rowSums(m), decreasing=TRUE)
word_data <- data.frame(word=names(word_freqs), freq=word_freqs)%>%
             mutate(word=ifelse(word=="marineprotectedarea","Marine Protected Area",word))%>%
             filter(!word %in% c("despite","can","also","used","using",
                                 "two","thus","showed","among","based","well","due","may","will",
                                 "less","three","found","use","especially","therefore","however","use","along",
                                 "often"))

# Generate Word Cloud

#saving this as a flattened image is a pain. I just did a screenshot of the R Studio preview
wordcloud2(word_data%>%filter(freq>10), size=1, color="random-dark", backgroundColor="white")

# Create a temporary directory
# temp_dir <- tempdir()
# temp_html <- file.path(temp_dir, "wordcloud.html")
# 
# saveWidget(wc, temp_html, selfcontained = FALSE)
# webshot(temp_html, "figures/eDNA_wordcloud.png", vwidth = 600, vheight = 600)
# 
# unlink(temp_dir, recursive = TRUE)

#round about way of saving - from ChatGPT = doesn't really work
# saveWidget(wc, "figures/wordcloud.html", selfcontained = FALSE)
# webshot("figures/wordcloud.html", "figures/eDNA_wordcloud.png", vwidth = 1000, vheight = 800)
