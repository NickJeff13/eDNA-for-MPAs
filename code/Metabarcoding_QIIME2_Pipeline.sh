#Processing eDNA metabarcode data using QIIME2 and cutadapt, followed by DADA2 taxonomy assignments using custom fish reference databases
#The qiime tutorials are useful and found at https://docs.qiime2.org/2022.2/tutorials/overview/#useful-points-for-beginners 
#First activate QIIME if it hasn't been, can also reactivate qiime if you close the window 
conda activate qiime2-2022.2

#check currently active conda environment
conda info
#activate tab completion
source tab-qiime
#Now import our data using a 'manifest' file of all fastq file names
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path pe33-16Smanifest \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path 16S-combined-demux.qza

#check out the data for visualization
qiime demux summarize \
  --i-data 16S-combined-demux.qza \
  --o-visualization 16S-demux-subsample.qzv
  
#Now trim primers
qiime cutadapt trim-paired \
--i-demultiplexed-sequences 16S-combined-demux.qza \
--p-front-f AGCGYAATCACTTGTCTYTTAA \
--p-front-r CRBGGTCGCCCCAACCRAA \
--p-error-rate 0.11 \
--output-dir  trimmed \
--verbose

#denoise using dada2 which infers ASVs 
#Note: using --p-n-threads = 0 will use all threads available 
qiime dada2 denoise-paired \
--i-demultiplexed-seqs trimmed/trimmed_sequences.qza \
--p-trim-left-f 10 \
--p-trim--left-r 10 \
--p-trunc-len-f  130 \
--p-trunc-len-r  130 \
--p-n-threads 0 \
--p-min-overlap 12 \
--p-pooling-method independent \
--output-dir trimmed/dada2out \
--verbose

#DADA2 didn't like the quality scores of my data (NovaSeq 6000 issue) so let's try merging reads with vsearch and denoising with DeBlur
qiime vsearch join-pairs \
--i-demultiplexed-seqs trimmed/trimmed_sequences.qza \
--o-joined-sequences joined_reads.qza \
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



