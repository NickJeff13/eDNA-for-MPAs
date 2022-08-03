#Obtain a custom "fish" 12S database using qiime and rescript
#txid77766 is Gnathostomata and includes sharks, rays, and all jawed fishes
qiime rescript get-ncbi-data \
    --p-query "txid7776[ORGN] AND (12S[Title] OR 12S ribosomal RNA[Title] OR 12S rRNA[Title]) AND (mitochondrion[Filter] OR plastid[Filter]) NOT environmental sample[Title] NOT environmental samples[Title] NOT environmental[Title] NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title]" \
    --p-ranks kingdom phylum class order family genus species \
    --p-rank-propagation \
    --p-n-jobs 20 \
    --o-sequences fish-12S-ref-seqs.qza \
    --o-taxonomy fish-12S-ref-tax.qza \
    --verbose
    
#Dereplicate reference data
qiime rescript dereplicate \
 --i-sequences fish-12S-ref-seqs.qza \
 --i-taxa fish-12S-ref-tax.qza \
 --p-mode 'uniq' \
 --p-threads 20 \
 --p-rank-handles 'disable' \
 --o-dereplicated-sequences fish-12S-ref-seqs-derep.qza \
 --o-dereplicated-taxa fish-12S-ref-tax-derep.qza
 
 #filter low-quality sequences
 qiime rescript cull-seqs \
 --i-sequences fish-12S-ref-seqs-derep.qza \
 --p-n-jobs 20 \
 --p-num-degenerates 5 \
 --p-homopolymer-length 8 \
 --o-clean-sequences fish-12S-ref-seqs-cull.qza
 
 #now filter by sequence length
 qiime rescript filter-seqs-length \
 --i-sequences fish-12S-ref-seqs-cull.qza \
 --p-global-min 90 \
 --p-global-max 1000 \
 --o-filtered-seqs fish-12S-ref-seqs-keep.qza \
 --o-discarded-seqs fish-12S-ref-seqs-discard.qza
