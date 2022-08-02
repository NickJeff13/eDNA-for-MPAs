###Running DADA2 on my metabarcoding data
#Modified from https://benjjneb.github.io/dada2/tutorial.html

#load libraries
library(dada2)


#these files are the output from cutadapt, after 12S primer removal
path<-"/hdd5/eDNA_Data/Raw/12S/Gully/output/"
list.files(path)


# Forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(path, pattern="_R1.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2.fastq", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 3)

#visualize quality profiles
plotQualityProfile(fnFs[1:2])
plotQualityProfile(fnRs[1:2])


# Place filtered files in filtered/ subdirectory
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                     truncLen=c(120,120),
                     maxN=0, maxEE=c(5,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)


errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
#plot error distributions
plotErrors(errF, nominalQ=TRUE)

#Now run the core DADA alogrithm
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)
#now merge forward and reverse reads to obtain full denoised sequences
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
# Inspect the merger data.frame from the first sample
head(mergers[[1]])


#Now construct our ASV table
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
table(nchar(getSequences(seqtab)))


#remove chimeras, which should be done after denoising to make identifying them easier
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)

#Now finally look at the number of reads that made it through each step of the pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

#Now let's assign taxonomy with our 12S database as the training set
taxa <- assignTaxonomy(seqtab.nochim, "../../../ReferenceData/12S_Actinopterygii_reflib_FEB25_2021.fa", multithread=TRUE)




