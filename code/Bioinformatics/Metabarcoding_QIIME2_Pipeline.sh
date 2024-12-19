#Processing eDNA metabarcode data using QIIME2 and cutadapt, followed by DADA2 taxonomy assignments using custom fish reference databases
#The qiime tutorials are useful and found at https://docs.qiime2.org/2022.2/tutorials/overview/#useful-points-for-beginners 
#First activate QIIME if it hasn't been, can also reactivate qiime if you close the window 
conda activate qiime2-2023.5 &&
source tab-qiime #activate tab completion

#check currently active conda environment
conda info

# View plugin, used to view any .qzv file in html and export tsv and fasta files
qiime tools view /path_to_file/filename.qzv

#load our params file which has location and marker/primer info
source /home/mcrg/Documents/Github/eDNA-for-MPAs/params

#Now import our data using a 'manifest' file of all fastq file names
#16S
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path pe33-16Smanifest \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path 16S-combined-demux.qza

#12S
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path pe33-12Smanifest \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path 12S-combined-demux.qza

#COI
qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path pe33-COImanifest \
--input-format PairedEndFastqManifestPhred33V2 \
--output-path COI-combined-demux.qza


#check out the data for visualization
#16S
qiime demux summarize \
  --i-data 16S-combined-demux.qza \
  --o-visualization 16S-demux-subsample.qzv ##save tsv file of per-sample-fastq-counts.tsv for optional step below ##
  
#12S
qiime demux summarize \
  --i-data 12S-combined-demux.qza \
  --o-visualization 12S-demux-subsample.qzv ##save tsv file of per-sample-fastq-counts.tsv for optional step below ##
    
  
 ## OPTIONAL: filter out samples with less than 100 reads (can set this to any number) ##
qiime demux filter-samples \
  --i-demux 16S-combined-demux.qza \
  --m-metadata-file /path_to_output_folder/per-sample-fastq-counts.tsv \
  --p-where 'CAST([forward sequence count] AS INT) > 100' \
  --o-filtered-demux /path_to_output_folder/filename_greater100reads.qza

#Now trim primers - COI 
qiime cutadapt trim-paired \
--i-demultiplexed-sequences COI-combined-demux.qza \
--p-cores 40 \
--p-front-f GGWACWGGWTGAACWGTWTAYCCYCC \
--p-front-r TAIACYTCIGGRTGICCRAARAAYCA \
--p-error-rate 0.11 \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-minimum-length 40 \
--o-trimmed-sequences COI-demux-trimmed.qza \
--output-dir  trimmed \
--verbose
#visualize the trimming results
qiime demux summarize --i-data COI-demux-trimmed.qza \
--o-visualization COI-trimmed-visual

#tTrim primers
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
--o-trimmed-sequences 16S-demux-trimmed.qza \
--output-dir  trimmed \
--verbose
#visualize the trimming results
qiime demux summarize --i-data 16s-demux-trimmed.qza \
--o-visualization 16s-trimmed-visual
':
=== Summary ===

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
':
MiFishU-F - these primers target ~180bp of the 12S rDNA region
5 GTCGGTAAAACTCGTGCCAGC 3
#MiFishU-R
3 GTTTGACCCTAATCTATGGGGTGATAC 5 > need to reverse this so its read 5 prime to 3 prime in cutadapt
' 
#12S -Mifish
qiime cutadapt trim-paired \
--i-demultiplexed-sequences 12S-combined-demux.qza \
--p-cores 40 \
--p-front-f GTCGGTAAAACTCGTGCCAGC \
--p-front-r NNNNNNCATAGTGGGGTATCTAATCCCAGTTTG \
--p-error-rate 0.15 \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-minimum-length 30 \
--o-trimmed-sequences 12s-demux-trimmed-2023-test2.qza \
--output-dir trimmed \
--verbose
&&
echo Trimming primer sequences complete!

#visualize the trimming results
qiime demux summarize --i-data 12s-demux-trimmed-2023-test2.qza \
--o-visualization 12s-trimmed-visual

#12S -248S-Forward with MiFishU-R
qiime cutadapt trim-paired \
--i-demultiplexed-sequences 12S-combined-demux.qza \
--p-cores 40 \
--p-front-f CGTGCCAGCCACCGCGGTT \
--p-front-r NNNNNNCATAGTGGGGTATCTAATCCCAGTTTG \
--p-error-rate 0.15 \
--p-discard-untrimmed \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--p-minimum-length 30 \
--o-trimmed-sequences 12s-Musq-demux-trimmed-2023-both.qza \
--output-dir trimmed \
--verbose
#visualize the trimming results
qiime demux summarize --i-data 12s-demux-trimmed-2023-test2.qza \
--o-visualization 12s-trimmed-visual

#denoise using dada2 which infers ASVs 
#Note: using --p-n-threads = 0 will use all threads available 
### for 16S use --p-trunc-len-f 125 and --p-trunc-len-r 125; 12S use 116 and 108 ###
# can add --p-min-overlap 12 or some other number if need be
#16S
qiime dada2 denoise-paired \
--i-demultiplexed-seqs 16S-demuxed-trimmed.qza \
--p-trunc-len-f  125 \
--p-trunc-len-r  125 \
--p-n-threads 0 \
--p-n-reads-learn 3000000 \
--p-pooling-method independent \
--output-dir dada2out-test \
--verbose

qiime dada2 denoise-paired \
--i-demultiplexed-seqs 16S-demuxed-trimmed.qza \
--p-trunc-len-f  145 \
--p-trunc-len-r  129 \
--p-n-threads 0 \
--p-n-reads-learn 3000000 \
--p-pooling-method independent \
--output-dir dada2out-test-2 \
--verbose
#trying longer with 16S as the sequences are 240bp long

qiime dada2 denoise-paired \
--i-demultiplexed-seqs 16S-combined-demux.qza \
--p-trunc-len-f  210 \
--p-trunc-len-r  210 \
--p-n-threads 0 \
--p-n-reads-learn 3000000 \
--p-pooling-method independent \
--output-dir dada2out \
--verbose

#12S - trying some different r-len truncs
qiime dada2 denoise-paired \
--i-demultiplexed-seqs 12S-demux-trimmed-2023.qza \
--p-trunc-len-f  189 \
--p-trunc-len-r  189 \
--p-n-threads 0 \
--p-min-overlap 8 \
--p-pooling-method independent \
--p-n-reads-learn 3000000 \
--output-dir denoised \
--verbose


#12S - trying some different r-len truncs
qiime dada2 denoise-single \
--i-demultiplexed-seqs 12s-demux-trimmed-2023.qza \
--p-trunc-len 125 \
--p-n-threads 0 \
--p-pooling-method independent \
--p-n-reads-learn 2000000 \
--output-dir ESIDenoisedSingle \
--verbose

#COI - trunc len 201 201 seems to work better than >220, and 191 191 did even better with some datasets
qiime dada2 denoise-paired \
--i-demultiplexed-seqs COI-demux-trimmed.qza \
--p-trunc-len-f  201 \
--p-trunc-len-r  201 \
--p-n-threads 0 \
--p-min-overlap 10 \
--p-pooling-method independent \
--p-n-reads-learn 3000000 \
--output-dir denoised \
--verbose

#Generate summaries of denoising stats and feature table

qiime feature-table summarize \
  --i-table denoised3/table.qza \
  --o-visualization denoised3/table.qzv \
  --m-sample-metadata-file ../../2023Perley-sample-metadata_ESI.tsv &&
qiime feature-table tabulate-seqs \
  --i-data denoisede3/representative_sequences.qza \
  --o-visualization denoised3/rep-seqs.qzv &&
qiime metadata tabulate \
  --m-input-file denoised/denoising_stats.qza \
  --o-visualization denoised/denoising-stats.qzv
  
#COI
qiime feature-table summarize \
  --i-table denoised3/table.qza \
  --o-visualization denoised3/table.qzv \
  --m-sample-metadata-file ../../seining2023-sample-metadata.tsv &&
qiime feature-table tabulate-seqs \
  --i-data dada2out/representative_sequences.qza \
  --o-visualization dada2out/rep-seqs.qzv &&
qiime metadata tabulate \
  --m-input-file dada2out/denoising_stats.qza \
  --o-visualization dada2out/denoising-stats.qzv
  
 qiime tools view /path_to_output_folder/filename_rep_seqs.qzv  ## export the ASV fasta file from the view for input into FuzzyID2 and BLAST

 #12S
 qiime feature-table summarize \
  --i-table ESIDenoised/table.qza \
  --o-visualization ESIDenoised/table.qzv \
  --m-sample-metadata-file ../2023Perley-sample-metadata_ESI.tsv &&
qiime feature-table tabulate-seqs \
  --i-data ESIDenoised/representative_sequences.qza \
  --o-visualization ESIDenoised/rep-seqs.qzv &&
qiime metadata tabulate \
  --m-input-file ESIDenoised/denoising_stats.qza \
  --o-visualization ESIDenoised/denoising-stats.qzv
  
   qiime feature-table summarize \
  --i-table ESIDenoisedSingle/table.qza \
  --o-visualization ESIDenoisedSingle/table.qzv \
  --m-sample-metadata-file ../2021-sample-metadata_ESIonly.tsv &&
qiime feature-table tabulate-seqs \
  --i-data ESIDenoisedSingle/representative_sequences.qza \
  --o-visualization ESIDenoisedSingle/rep-seqs.qzv &&
qiime metadata tabulate \
  --m-input-file ESIDenoisedSingle/denoising_stats.qza \
  --o-visualization ESIDenoisedSingle/denoising-stats.qzv
  
  qiime tools view /path_to_output_folder/filename_rep_seqs.qzv  ## export the ASV fasta file from the view for input into FuzzyID2 and BLAST

 ### export results to biom formatted file
qiime tools export \
--input-path dada2out-test/table.qza \
--output-path dada2out-test/ESI16S_filtered_table_biom ##specifying a folder output here, this tool will automatically export a file called 'feature-table.biom' to this folder

### convert biom to tsv
biom convert -i dada2out-test/ESI16S_filtered_table_biom/feature-table.biom \
-o dada2out-test/ESI16S_filtered_table_biom/ESI16S_feature_table_export.tsv \
--to-tsv

### OPTIONAL filtering after exporting to tsv
## Remove rare ASV's by calculating if an ASV has a read number that is less than 0.1% of the total read number of that ASV across all samples. 
## This is summing across columns in the exported feature table, calculating 0.1% of that sum, and removing all instances where read numbers were less than that number.
 
 #Generate a phylogenetic tree from our data
 cd dada2out/
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
  --m-metadata-file ../../../2023Perley-sample-metadata_ESI.tsv \
  --output-dir COI-core-metrics-results
 
####################
######TAXONOMY######
#####################
## there are lots of ways to do taxonomy, including just blasting, or building a reference database using rescript (below), or using FuzzyID2 with a custom library
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
  
:'OLD CODE and ISSUES: 
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


