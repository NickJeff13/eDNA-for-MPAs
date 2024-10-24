#Processing eDNA metabarcode data using QIIME2 and cutadapt, followed by DADA2 taxonomy assignments using custom fish reference databases
#The qiime tutorials are useful and found at https://docs.qiime2.org/2022.2/tutorials/overview/#useful-points-for-beginners 
#First activate QIIME if it hasn't been, can also reactivate qiime if you close the window 
conda activate qiime2_2023_5 &&
source tab-qiime #activate tab completion

#check currently active conda environment
conda info

# View plugin, used to view any .qzv file in html and export tsv and fasta files
qiime tools view /path_to_file/filename.qzv

####################
#########12S########
####################

#####################
#STEP 1 - Import Data#
#####################

#create a manifest file maually with "sample-id","forward-absolute-filepath","reverse-absolute-filepath", save without a file extension but as a tsv
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path pe33-12Smanifest-Musquash \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path Musquash-12S-combined-demux.qza  #input-path is the manifest .tsv file listing the full path to each .gz file

#check out the data for visualization
qiime demux summarize \
  --i-data /mnt/sda2/eDNA/Musquash/Data/12S/Musquash-12S-combined-demux.qza \
  --o-visualization /mnt/sda2/eDNA/Musquash/Data/12S/Musquash-12S-combined-demux.qzv ##save tsv file of per-sample-fastq-counts.tsv for optional step below ##

## OPTIONAL: filter out samples with less than 100 reads (can set this to any number) ##
qiime demux filter-samples \
  --i-demux 16S-combined-demux.qza \
  --m-metadata-file /path_to_output_folder/per-sample-fastq-counts.tsv \
  --p-where 'CAST([forward sequence count] AS INT) > 100' \
  --o-filtered-demux /path_to_output_folder/filename_greater100reads.qza

#####################  
#STEP 2 - Trim Primers#
#####################
#See FishPrimers.txt for list of primers to use  
qiime cutadapt trim-paired \
--i-demultiplexed-sequences Musquash-12S-combined-demux.qza \
--p-cores 60 \
--p-front-f NNNNNNGTCGGTAAAACTCGTGCCAGC \
--p-front-r NNNNNNCATAGTGGGGTATCTAATCCCAGTTTG \
--p-error-rate 0.15 \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-minimum-length 30 \
--o-trimmed-sequences Musquash-12S-combined-demux-trimmed.qza \
--output-dir trimmed \
--verbose
#visualize the trimming results
qiime demux summarize --i-data Musquash-12S-combined-demux-trimmed.qza \
--o-visualization Musquash-12S-combined-demux-trimmed.qzv 

#####################
##STEP 3 - Denoise using dada2##
#####################
#denoise using dada2 which infers ASVs
#Note: using --p-n-threads = 0 will use all threads available 
### for 16S use --p-trunc-len-f 125 and --p-trunc-len-r 125
#for 12S use --p-trunc-len-f 116 and --p-trunc-len-r 108
#for COI use --p-trunc-len-f 180 and --p-trunc-len-r 180 ###
# can add --p-min-overlap 12 or some other number if need be

qiime dada2 denoise-paired \
--i-demultiplexed-seqs Musquash-12S-combined-demux-trimmed.qza \
--p-trunc-len-f  116 \
--p-trunc-len-r  108 \
--p-n-threads 0 \
--p-pooling-method independent \
--output-dir dada2out_12S \
--verbose

#arguments not used
#how many reads to use in training the algorithm
--p-n-reads-learn 3000000 \ 
#merge overlap
--p-min-overlap 8 \

#####################
##STEP 4 - Generate summaries##
#####################
#Generate summaries of denoising stats and feature table
qiime feature-table summarize \
  --i-table dada2out_12S/table.qza \
  --o-visualization dada2out_12S/table.qzv \
  --m-sample-metadata-file Musquash-12S-metadata_dada2.tsv &&
qiime feature-table tabulate-seqs \
  --i-data dada2out_12S/representative_sequences.qza \
  --o-visualization dada2out_12S/representative_sequences.qzv &&
qiime metadata tabulate \
  --m-input-file dada2out_12S/denoising_stats.qza \
  --o-visualization dada2out_12S/denoising_stats.qzv
  
qiime tools view dada2out_12S/table.qzv 
qiime tools view dada2out_12S/representative_sequences.qzv  ## export the ASV fasta file from the representative sequences view for input into FuzzyID2 and BLAST
qiime tools view dada2out_12S/denoising_stats.qzv  ## export the table to compare read loss through filtering steps


#####################
##STEP 5 - Visualize Results##
#####################
#export results to biom formatted file
qiime tools export \
--input-path dada2out_12S/table.qza \
--output-path dada2out_12S/Musquash_12S_filtered_table_biom ##specifying a folder output here, this tool will automatically export a file called 'feature-table.biom' to this folder

#convert biom to tsv
biom convert -i dada2out_12S/Musquash_12S_filtered_table_biom/feature-table.biom \
-o dada2out_12S/Musquash_12S_filtered_table_biom/Musquash_12S_feature_table_export.tsv \
--to-tsv

#Generate a phylogenetic tree from our data
cd dada2out_12S/
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences representative_sequences.qza \
--o-alignment aligned-rep-seqs.qza \
--o-masked-alignment masked-aligned-rep-seqs.qza \
--o-tree unrooted-tree.qza \
--o-rooted-tree rooted-tree.qza

#now use the rooted tree to generate some biodiversity stats
qiime diversity core-metrics-phylogenetic \
--i-phylogeny rooted-tree.qza \
--i-table table.qza \
--p-sampling-depth 1500 \
--p-n-jobs-or-threads auto \
--m-metadata-file ../Musquash-12S-metadata_dada2.tsv \
--output-dir 12S-core-metrics-results

#####################
##STEP 6 - Taxonomy##
#####################
#Next, calculate the taxonomy for the ASVs using the code below, but much of the filtering and consolidating will be done in R
#consolidate the taxonomy for the ASVs:
  #manually BLAST anything that's unassigned
  #assign things to the lowest common taxon if multiple mataches

## there are lots of ways to do taxonomy, including just blasting, or building a reference database using rescript (see rescript_createReferenceDB.sh), or using FuzzyID2 with a custom library

###BLAST###
#using the custom 12S database built with rescript --> rescript_createReferenceDB.sh
#for 12S, anything under ~97% perc. ident may not be valid down to species level
qiime feature-classifier blast \
--i-query representative_sequences.qza \
--i-reference-reads ../../../DBs/fish-12S-ref-seqs-FINAL.qza \
--p-maxaccepts 20 \
--p-perc-identity 0.9 \
--o-search-results Musquash_12S_blastOutput.qza \
--o-visualization Musquash_12S_blastOutput.qzv

qiime feature-classifier classify-consensus-blast \
--i-query representative_sequences.qza \
--i-reference-reads ../../../DBs/fish-12S-ref-seqs-FINAL.qza \
--i-reference-taxonomy ../../../DBs/fish-12S-ref-tax-FINAL.qza \
--p-maxaccepts 20 \
--p-perc-identity 0.9 \
--o-search-results Musquash_12S_blastclassifierOutput.qza \
--o-classification Musquash_12S_blastclassifierOutput_tax.qza \
--output-dir BlastClassifierTaxonomy


#once taxonomy is calculated and finalized, use R to calculate biodiversity stats, etc
##see eDNA_BiodiversityStats_AfterQIIME2.R

#####################
#########COI#########
#####################

#####################
##STEP 1 - Import Data##
#####################

#create a manifest file maually with "sample-id","forward-absolute-filepath","reverse-absolute-filepath", save without a file extension but as a tsv
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path pe33-COImanifest-Musquash \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path Musquash-COI-combined-demux.qza

#check out the data for visualization
qiime demux summarize \
  --i-data Musquash-COI-combined-demux.qza \
  --o-visualization Musquash-COI-combined-demux.qzv ##save tsv file of per-sample-fastq-counts.tsv for optional step below ##
    
## OPTIONAL: filter out samples with less than 100 reads (can set this to any number) ##
qiime demux filter-samples \
  --i-demux 16S-combined-demux.qza \
  --m-metadata-file /path_to_output_folder/per-sample-fastq-counts.tsv \
  --p-where 'CAST([forward sequence count] AS INT) > 100' \
  --o-filtered-demux /path_to_output_folder/filename_greater100reads.qza

#####################
##STEP 2 - Trim Primers##
#####################
#See FishPrimers.txt for list of primers to use
qiime cutadapt trim-paired \
--i-demultiplexed-sequences Musquash-COI-combined-demux.qza \
--p-cores 60 \
--p-front-f GGWACWGGWTGAACWGTWTAYCCYCC \
--p-front-r GGRGGRTASACSGTTCASCCSGTSCC \
--p-error-rate 0.15 \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-minimum-length 30 \
--o-trimmed-sequences Musquash-COI-combined-demux-trimmed.qza \
--output-dir trimmed \
--verbose
#visualize the trimming results
qiime demux summarize --i-data Musquash-COI-combined-demux-trimmed.qza \
--o-visualization Musquash-COI-combined-demux-trimmed.qzv

qiime tools view /mnt/sda2/eDNA/Musquash/Data/COI/Musquash-COI-combined-demux-trimmed.qzv

#####################
##STEP 3 - Denoise using dada2##
#####################
#denoise using dada2 which infers ASVs
#Note: using --p-n-threads = 0 will use all threads available 
### for 16S use --p-trunc-len-f 125 and --p-trunc-len-r 125
#for 12S use --p-trunc-len-f 116 and --p-trunc-len-r 108
#for COI use --p-trunc-len-f 180 and --p-trunc-len-r 180 ###
# can add --p-min-overlap 12 or some other number if need be

qiime dada2 denoise-paired \
--i-demultiplexed-seqs Musquash-COI-combined-demux-trimmed.qza \
--p-trunc-len-f  180 \
--p-trunc-len-r  180 \
--p-n-threads 0 \
--p-pooling-method independent \
--output-dir dada2out_COI \
--verbose


#####################
##STEP 4 - Generate summaries##
#####################
#Generate summaries of denoising stats and feature table

qiime feature-table summarize \
  --i-table /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/table.qza \
  --o-visualization /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/table.qzv \
  --m-sample-metadata-file /mnt/sda2/eDNA/Musquash/Data/COI/Musquash-COI-metadata_dada2.tsv &&
qiime feature-table tabulate-seqs \
  --i-data /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/representative_sequences.qza \
  --o-visualization /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/representative_sequences.qzv &&
qiime metadata tabulate \
  --m-input-file /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/denoising_stats.qza \
  --o-visualization /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/denoising_stats.qzv
  
qiime tools view /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/table.qzv 
qiime tools view /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/representative_sequences.qzv  ## export the ASV fasta file from the representative sequences view for input into BLAST
qiime tools view /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/denoising_stats.qzv  ## export the table to compare read loss through filtering steps

#####################
##STEP 5 - Visualize Results##
#####################
#export results to biom formatted file#
qiime tools export \
--input-path /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/table.qza \
--output-path /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/Musquash_COI_filtered_table_biom ##specifying a folder output here, this tool will automatically export a file called 'feature-table.biom' to this folder

#convert biom to tsv
biom convert -i /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/Musquash_COI_filtered_table_biom/feature-table.biom \
-o /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/Musquash_COI_filtered_table_biom/Musquash_COI_feature_table_export.tsv \
--to-tsv

#Generate a phylogenetic tree from our data
qiime phylogeny align-to-tree-mafft-fasttree \
--i-sequences /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/representative_sequences.qza \
--o-alignment /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/aligned-rep-seqs.qza \
--o-masked-alignment /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/masked-aligned-rep-seqs.qza \
--o-tree /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/unrooted-tree.qza \
--o-rooted-tree /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/rooted-tree.qza

#now use the rooted tree to generate some biodiversity stats
qiime diversity core-metrics-phylogenetic \
--i-phylogeny /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/rooted-tree.qza \
--i-table /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/table.qza \
--p-sampling-depth 1500 \
--p-n-jobs-or-threads auto \
--m-metadata-file /mnt/sda2/eDNA/Musquash/Data/COI/Musquash-COI-metadata_dada2.tsv \
--output-dir /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/COI-core-metrics-results

#####################
##STEP 6 - Taxonomy##
#####################
#Next, calculate the taxonomy for the ASVs using the code below, but much of the filtering and consolidating will be done in R
#consolidate the taxonomy for the ASVs:
  #manually BLAST anything that's unassigned
  #assign things to the lowest common taxon if multiple mataches

#once taxonomy is calculated and finalized, use R to calculate biodiversity stats, etc
#see eDNA_BiodiversityStats_AfterQIIME2.R

# there are lots of ways to do taxonomy, including just blasting, or building a reference database using rescript (see rescript_createReferenceDB.sh), or using FuzzyID2 with a custom library
 
###BLAST###
#using the custom COI database built with rescript --> rescript_createReferenceDB.sh
#for COI, anything under ~97% perc. ident may not be valid down to species level
#blast only
qiime feature-classifier blast \
--i-query representative_sequences.qza \
--i-reference-reads /home/ABL/eDNA/Musquash/DBs/COI/fish-COI-ref-seqs-19Sept2023-FINAL.qza \
--p-maxaccepts 20 \
--p-perc-identity 0.9 \
--o-search-results BlastTaxonomy_DB19Sept2023/Musquash_COI_blastOutput-19Sept2023.qza \
--verbose

#export results
qiime tools export \
--input-path /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/BlastTaxonomy_DB19Sept2023/Musquash_COI_blastOutput-19Sept2023.qza \
--output-path /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/BlastTaxonomy_DB19Sept2023/Musquash_COI_blastOutput-19Sept2023 ##specifying a folder output here, 

# #blast and taxonomy consensus - not really useful
# qiime feature-classifier classify-consensus-blast \
# --i-query representative_sequences.qza \
# --i-reference-reads /home/ABL/eDNA/Musquash/DBs/COI/fish-COI-ref-seqs-25Aug2023-FINAL.qza \
# --i-reference-taxonomy /home/ABL/eDNA/Musquash/DBs/COI/fish-COI-ref-tax-25Aug2023-FINAL.qza \
# --p-maxaccepts 10 \
# --p-perc-identity 0.7 \
# --o-search-results /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/BlastClassifierTaxonomy_DB25Aug2023/Musquash_COI_blastclassifierOutput-19Sept2023.qza \
# --o-classification /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/BlastClassifierTaxonomy_DB25Aug2023/Musquash_COI_blastclassifierOutput_tax-19Sept2023.qza
# 
# 
# #export results
# qiime tools export \
# --input-path /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/BlastClassifierTaxonomy_DB25Aug2023/Musquash_COI_blastclassifierOutput-19Sept2023.qza \
# --output-path /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/BlastClassifierTaxonomy_DB25Aug2023/Musquash_COI_blastclassifierOutput-19Sept2023 ##specifying a folder output here, 
# 
# qiime tools export \
# --input-path /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/BlastClassifierTaxonomy_DB25Aug2023/Musquash_COI_blastclassifierOutput_tax-19Sept2023.qza \
# --output-path /mnt/sda2/eDNA/Musquash/Data/COI/dada2out_COI/BlastClassifierTaxonomy_DB25Aug2023/Musquash_COI_blastclassifierOutput_tax-19Sept2023 ##specifying a folder output here, 
# 

############################################
##################EXTRA CODE################
############################################

#####Now trim primers#####
#See FishPrimers.txt for list of primers to use
#16S#
qiime cutadapt trim-paired \
--i-demultiplexed-sequences 16S-combined-demux.qza \
--p-cores 40 \
--p-front-f AGCGYAATCACTTGTCTYTTAA \
--p-front-r CRBGGTCGCCCCAACCRAA \
--p-error-rate 0.11 \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-minimum-length 40 \
--o-trimmed-sequences 16sGul-demux-trimmed.qza \
--output-dir  trimmed \
--verbose
#visualize the trimming results
qiime demux summarize --i-data 16s-demux-trimmed.qza \
--o-visualization 16s-trimmed-visual
#Example output:
'=== Summary ===

Total read pairs processed:          5,256,200
  Read 1 with adapter:               3,544,875 (67.4%)
  Read 2 with adapter:               3,531,749 (67.2%)

== Read fate breakdown ==
Pairs that were too short:                   0 (0.0%)
Pairs discarded as untrimmed:        1,812,812 (34.5%)
Pairs written (passing filters):     3,443,388 (65.5%)

Total basepairs processed: 1,650,284,162 bp
  Read 1:   825,252,162 bp
  Read 2:   825,032,000 bp
Total written (filtered):    943,244,238 bp (57.2%)
  Read 1:   466,537,198 bp
  Read 2:   476,707,040 bp
'





#Next, calculate the taxonomy for the ASVs using the code below, but much of the filtering and consolidating will be done in R
#consolidate the taxonomy for the ASVs:
  #manually BLAST anything that's unassigned
  #assign things to the lowest common taxon if multiple mataches


#once taxonomy is calculated and finalized, use R to calculate biodiversity stats, etc
###Necessary filtering:
## Remove rare ASV's by calculating if an ASV has a read number that is less than 1% of the total read number of that ASV across all samples. 
    ## This is summing across columns in the exported feature table, calculating 1% of that sum, and removing all instances where read numbers were less than that number.
## Remove ASVs that contain <1% of the total reads of that entire sample across all ASVs
## This is summing across columns in the exported feature table, calculating 1% of that sum, and removing all instances where read numbers were less than that number.
##see eDNA_BiodiversityStats_AfterQIIME2.R


 
####################
######TAXONOMY######
#####################
## there are lots of ways to do taxonomy, including just blasting, or building a reference database using rescript (see rescript_createReferenceDB.sh), or using FuzzyID2 with a custom library
#####rescript#####
#using rescript to train our classifier
 #Build and evaluate classifier
 #here, the --o-classifier output is of type TaxonomicClassifier and the -o-observed-taxonomy is FeatureData[Taxonomy] (same as --i-taxonomy)
 qiime rescript evaluate-fit-classifier \
 --i-sequences fish-12S-ref-seqs-FINAL.qza \
 --i-taxonomy fish-12S-ref-tax-FINAL.qza \
 --p-n-jobs -1 \
 --o-classifier ncbi-12S-fish-refseqs-classifier.qza \
 --o-evaluation ncbi-12S-fish-refseqs-classifier-evaluation.qzv \
 --o-observed-taxonomy ncbi-12S-fish-refseqs-predicted-taxonomy.qza \
 --output-dir 12S-Classifier \
 --verbose 

 # FutureWarning: Passing a set as an indexer is deprecated and will raise in a future version. Use a list instead.
 #  taxa = taxa.loc[seq_ids]
 #UserWarning: The TaxonomicClassifier artifact that results from this method was trained using scikit-learn version 0.24.1. It cannot be used with other versions of scikit-learn. (While the classifier may complete successfully, the results will be unreliable.)
 
 
 qiime rescript evaluate-taxonomy \
 --i-taxonomies fish-16S-ref-taxa-keep.qza ncbi-16S-fish-refseqs-predicted-taxonomy.qza \
 --p-labels ref-taxonomy predicted-taxonomy \
 --o-taxonomy-stats 16S-ref-taxonomy-evaluation.qzv \
 --verbose
 
#Now back to qiime to do our taxonomy
  qiime feature-classifier classify-sklearn \
  --i-classifier ../../../../ReferenceData/ncbi-16S-fish-refseqs-classifier.qza \
  --i-reads representative_sequences.qza \
  --o-classification 16S-taxonomy.qza

qiime metadata tabulate \
  --m-input-file 16S-taxonomy.qza \
  --o-visualization 16S-taxonomy.qzv

#######16S############


###rescript###

#using rescript to train our classifier
qiime rescript filter-taxa \
--i-taxonomy fish-16S-ref-tax.qza \
--m-ids-to-keep-file fish-12S-ref-seqs-FINAL.qza \
--o-filtered-taxonomy fish-16S-ref-taxa-FINAL.qza

qiime rescript evaluate-taxonomy \
--i-taxonomies fish-16S-ref-taxa-keep.qza \
--o-taxonomy-stats fish-16S-ref-tax-keep-eval.qzv

qiime metadata tabulate \
--m-input-file fish-16S-ref-taxa-keep.qza \
--o-visualization fish-16S-ref-tax-keep.qzv &&
qiime rescript evaluate-seqs \
--i-sequences fish-16S-ref-seqs-keep.qza \
--p-kmer-lengths 32 16 8 \
--o-visualization fish-16S-ref-seqs-keep-eval.qzv
 
 #Build and evaluate classifier
 #here, the --o-classifier output is of type TaxonomicClassifier and the -o-observed-taxonomy is FeatureData[Taxonomy] (same as --i-taxonomy)
 qiime rescript evaluate-fit-classifier \
 --i-sequences fish-16S-ref-seqs-keep.qza \
 --i-taxonomy fish-16S-ref-taxa-keep.qza \
 --p-n-jobs -1 \
 --o-classifier ncbi-16S-fish-refseqs-classifier.qza \
 --o-evaluation ncbi-16S-fish-refseqs-classifier-evaluation.qzv \
 --o-observed-taxonomy ncbi-16S-fish-refseqs-predicted-taxonomy.qza \
 --output-dir 16S-Classifier \
 --verbose 
 
 qiime rescript evaluate-taxonomy \
 --i-taxonomies fish-16S-ref-taxa-keep.qza ncbi-16S-fish-refseqs-predicted-taxonomy.qza \
 --p-labels ref-taxonomy predicted-taxonomy \
 --o-taxonomy-stats 16S-ref-taxonomy-evaluation.qzv \
 --verbose
 
#Now back to qiime to do our taxonomy
  qiime feature-classifier classify-sklearn \
  --i-classifier ../../../../ReferenceData/ncbi-16S-fish-refseqs-classifier.qza \
  --i-reads representative_sequences.qza \
  --o-classification 16S-taxonomy.qza

qiime metadata tabulate \
  --m-input-file 16S-taxonomy.qza \
  --o-visualization 16S-taxonomy.qzv
  
#This gave us a long tsv file for download which shows the features and assigned taxonomy with confidence values
#Now we'll use ANCOM to test differences among sites
#First subset by just the BrownsBank samples
qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file ../../../2021-sample-metadata.tsv \
  --p-where "[location]='BrownsBank'" \
  --o-filtered-table BBL-table.qza
  
  qiime composition add-pseudocount \
  --i-table BBL-table.qza \
  --o-composition-table comp-BBL-table.qza
  
  #this can only be done comparing two groups - can't have all the same or more than 2 unique locations,sites, etc
  qiime composition ancom \
  --i-table comp-BBL-table.qza \
  --m-metadata-file ../../../2021-sample-metadata.tsv \
  --m-metadata-column location \
  --o-visualization 16SBBL-ancom-subject.qzv
  
:OLD CODE and ISSUES: 
DADA2 didnt like the quality scores of my data (NovaSeq 6000 issue) so lets try merging reads with vsearch and denoising with DeBlur
Actually not entirely true - it runs, though error plots look weird because of NovaSeq quality score binning. Regardless, shortening my sequences to 130bp seems to have worked for dada2.
Below are commands to run vsearch to join pairs and then deblur 
qiime vsearch join-pairs \
--i-demultiplexed-seqs 12S-demux-trimmed.qza \
--p-minlen 100 \
--o-joined-sequences joined-trimmed-reads.qza \
--verbose \
--p-threads 8 #8 threads is the max here

#Check out quality scores
qiime quality-filter q-score \
      --i-demux joined_reads.qza \
      --o-filtered-sequences joined-filtered.qza \
      --o-filter-stats joined-filter-stats.qza

qiime deblur denoise-16S \
      --i-demultiplexed-seqs joined-filtered.qza \
      --p-trim-length 160 \
      --p-sample-stats \
      --p-jobs-to-start 32 \
      --o-stats deblur-stats.qza \
      --o-representative-sequences rep-seqs-deblur.qza \
      --o-table table-deblur.qza \
      --verbose
      
#Visualize deblurring
qiime deblur visualize-stats \
       --i-deblur-stats deblur-stats.qza \
       --o-visualization deblur-stats.qzv

qiime feature-table tabulate-seqs \
       --i-data rep-seqs-deblur.qza \
       --o-visualization rep-seqs-deblur.qzv

qiime feature-table summarize \
       --i-table table-deblur.qza \
       --m-sample-metadata-file metadata.tsv \
       --o-visualization table-deblur.qzv
'


