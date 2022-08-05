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

#Generate summaries of denoising stats and feature table
qiime feature-table summarize \
  --i-table trimmed/dada2out/table.qza \
  --o-visualization trimmed/dada2out/table.qzv \
  --m-sample-metadata-file ../2021-sample-metadata.tsv &&
qiime feature-table tabulate-seqs \
  --i-data trimmed/dada2out/representative_sequences.qza \
  --o-visualization trimmed/dada2out/rep-seqs.qzv &&
qiime metadata tabulate \
  --m-input-file trimmed/dada2out/denoising-stats.qza \
  --o-visualization trimmed/dada2out/denoising-stats.qzv
 
 #Generate a phylogenetic tree from our data
 cd trimmed/dada2out/
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
  --m-metadata-file ../../../2021-sample-metadata.tsv \
  --output-dir core-metrics-results
  
  ###TAXONOMY####
  #using rescript to train our classifier
  qiime rescript filter-taxa \
  --i-taxonomy fish-16S-ref-tax.qza \
  --m-ids-to-keep-file fish-16S-ref-seqs-keep.qza \
  --o-filtered-taxonomy fish-16S-ref-taxa-keep.qza
  
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
  
  
  
  
:' DADA2 didnt like the quality scores of my data (NovaSeq 6000 issue) so lets try merging reads with vsearch and denoising with DeBlur
Actually not entirely true - it runs, though error plots look weird because of NovaSeq quality score binning. Regardless, shortening my sequences to 130bp seems to have worked for dada2.
Below are commands to run vsearch to join pairs and then deblur 
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
'


