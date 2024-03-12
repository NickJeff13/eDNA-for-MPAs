#Loading required packages for the phylogenetic analysis
#libraries - install first! 
library(seqinr)
library(ape)
library(msa) # this package is available through the Bioconductor packages. If you do not have the BiocManager package installed, you could install this using the following lines of code... then load it using library(msa)

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")

#BiocManager::install("msa")
#library(msa)

### Note: this becomes quite a bit of data fast - almost 1Gb of alignment objects
## Read sequences from FASTA file of the ASVs in the output folder
## A new note: this tree looked HORRIBLE because there's so many ASVs, so I am going to try to subset the fasta file with the ones we used 

asvs <- read.table("data/2022Data/SAB/COI/COI_feature_table_export.tsv", header = T)
asvs.keep <- asvs$OTU.ID #gives us just the ASV names from dada2

sequence <- readDNAStringSet("data/2022Data/SAB/COI/dna-sequences.fasta")

sequence.filt <- sequence[1:500]

## Perform multiple sequence alignment using the ClustalW algorithm

my_alignment <- msa(sequence.filt, "ClustalW")

## Visualize the alignment if you want - doesn't seem to work for Nick but you can view it just be typing out my_alignment
#msaPrettyPrint(my_alignment, output="dvi",
 #              showNames="left", showLogo="none", askForOverwrite=FALSE)

## Compute distance matrix
my_alignment_sequence <- msaConvert(my_alignment, type="seqinr::alignment")

distance_alignment <- dist.alignment(my_alignment_sequence)

## compute phylogenetic tree using neighbor joining

Tree <- bionj(distance_alignment)
#replace the tree tip labels with species names
Tree$tip.label <- asvs$V26

## display phylogenetic tree - in this case 'plot' is short for plot.phylo but we can just say plot because "Tree" is class 'phylo', so R knows how to plot it

plot(Tree, type="phylogram", align.tip.label=T)
plot(Tree, type="fan", show.tip.label=F, edge.width=2, cex=0.7, no.margin=F)
plot(Tree, type="unrooted", use.edge.length=T, show.tip.label=T, no.margin=T)

save.image("data/Tree.RData")
