setwd("~/eDNA/Musquash/Data/12S/dada2out_12S_Test4")

library("ape")
library("ggplot2")

# Bioconductor version
library(BiocManager)
BiocManager::install("ggtree")
BiocManager::install("Biostrings")

#load metadata
metadata <- read.delim("~/eDNA/Musquash/Data/12S/Musquash-12S-metadata_dada2.tsv", sep="\t")

#load tree file
nwk <- system.file("extdata", "Musquash12S_UnrootedTree_QIIME2.nwk", package="ggtree")
tree <- read.tree("Musquash12S_UnrootedTree_QIIME2.nwk")
ggtree(tree)
