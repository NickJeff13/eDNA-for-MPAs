#Processing ARMS scrapings/bulk sample metabarcode data using QIIME2 and cutadapt, followed by DADA2 taxonomy assignments using custom fish reference databases
#The qiime tutorials are useful and found at https://docs.qiime2.org/2022.2/tutorials/overview/#useful-points-for-beginners 

#First activate QIIME if it hasn't been, can also reactivate qiime if you close the window 

conda activate qiime2-amplicon-2024.10 &&
source tab-qiime #activate tab completion


#Now import our data using a 'manifest' file of all fastq file names
#COI
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path pe33-ARMSmanifest \
--input-format CasavaOneEightSingleLanePerSampleDirFmt \
--output-path ARMS-combined-demux.qza


#check out the data for visualization
qiime demux summarize \
  --i-data ARMS-combined-demux.qza \
  --o-visualization ARMS-demux-subsample.qzv ##save tsv file of per-sample-fastq-counts.tsv for optional step below ##
  
  
  
#Now trim primers - COI Leray primers
  
qiime cutadapt trim-paired \
--i-demultiplexed-sequences ARMS-combined-demux.qza \
--p-cores 40 \
--p-front-f GGWACWGGWTGAACWGTWTAYCCYCC \
--p-front-r TAIACYTCIGGRTGICCRAARAAYCA \
--p-error-rate 0.11 \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-minimum-length 40 \
--o-trimmed-sequences ARMS-demux-trimmed.qza \
--verbose

#visualize the trimming results
qiime demux summarize --i-data ARMS-demux-trimmed.qza \
--o-visualization ARMS-trimmed-visual

#Denoise using DADA2

qiime dada2 denoise-paired \
--i-demultiplexed-seqs ARMS-demux-trimmed.qza \
--p-trunc-len-f  230 \
--p-trunc-len-r  230 \
--p-n-threads 0 \
--p-min-overlap 12 \
--p-max-ee-f 5 \
--p-max-ee-r 5 \
--p-pooling-method independent \
--p-n-reads-learn 2000000 \
--output-dir denoised \
--verbose

#Check denoising stats
qiime metadata tabulate \
  --m-input-file denoised/denoising_stats.qza \
  --o-visualization denoised/denoising-stats.qzv
  
#Next generate ASV table
qiime feature-table summarize \
  --i-table denoised/table.qza \
  --o-visualization denoised/table.qzv \
  --m-sample-metadata-file ../ARMS_metabarcoding_metadata.tsv &&
qiime feature-table tabulate-seqs \
  --i-data denoised/representative_sequences.qza \
  --o-visualization denoised/rep-seqs.qzv 
  
  
  
  
 ### export results to biom formatted file
qiime tools export \
--input-path denoised/table.qza \
--output-path denoised/ARMS_metabarcoding_filtered_table_biom ##specifying a folder output here, this tool will automatically export a file called 'feature-table.biom' to this folder

### convert biom to tsv
biom convert -i denoised/ARMS_metabarcoding_filtered_table_biom/feature-table.biom \
-o denoised/ARMS_metabarcoding_filtered_table_biom/ARMS_metabar_feature_table_export.tsv \
--to-tsv

 #Generate a phylogenetic tree from our data
 cd denoised/
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
  --m-metadata-file ../../ARMS_metabarcoding_metadata.tsv \
  --output-dir core-metrics-results
