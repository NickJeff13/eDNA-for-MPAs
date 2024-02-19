###Running DADA2 on my metabarcoding data
#Modified from https://benjjneb.github.io/dada2/bigdata.html

#load libraries
library(dada2)


#these files are the output from cutadapt, after 12S primer removal
path<-"/hdd5/eDNA_Data/Raw/12S/cutadapt"
list.files(path)


# Try using only forward reads for 12S - many of my pairs don't overlap properly
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 2) #the three here means the third underscore in the filename 

#visualize quality profiles
plotQualityProfile(fnFs[1:2])

plotComplexity(fnFs[1:2])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
names(filtFs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, 
                     truncLen=126,
                     maxN=0, maxEE=4, truncQ=10, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)


errF <- learnErrors(filtFs, nbases=1e+10, multithread=TRUE, verbose=1)
#plot error distributions
plotErrors(errF, nominalQ=TRUE)

#Now run the core DADA algorithm
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)



#Now construct our ASV table
seqtab <- makeSequenceTable(dadaFs)
dim(seqtab)
table(nchar(getSequences(seqtab)))



#Now finally look at the number of reads that made it through each step of the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF")
rownames(track) <- sample.names
head(track)

#Now let's assign taxonomy with our 12S database as the training set
#take a look at an example reference for assignSpecies

twelveS.seqs <- getSequences(object = seqtab)

