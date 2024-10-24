############################################
#############12S############################
############################################
#Obtain a custom "fish" 12S database using qiime and rescript
#txid77766 is Gnathostomata and includes sharks, rays, and all jawed fishes
#activate the rescript environment
conda activate rescript

qiime rescript get-ncbi-data \
    --p-query "txid7776[ORGN] AND\
    (12S[Title] OR 12S ribosomal RNA[Title] OR 12S rRNA[Title] OR 12S mitochondrial[Title] OR mitochondrion[Title] OR small subunit ribosomal[Title]) AND\
    (mitochondrion[Filter] OR plastid[Filter] OR genomic dna[Filter]) NOT\
    environmental sample[Title] NOT\
    environmental samples[Title] NOT\
    environmental[Title] NOT\
    uncultured[Title] NOT\
    unclassified[Title] NOT\
    unidentified[Title] NOT\
    unverified[Title]" \
    --p-ranks kingdom phylum class order family genus species \
    --p-rank-propagation \
    --p-n-jobs 20 \
    --o-sequences fish-12S-ref-seqs_24Aug2023.qza \
    --o-taxonomy fish-12S-ref-tax_24Aug2023.qza \
    --verbose
    
#Dereplicate reference data - sequences and taxonomy
qiime rescript dereplicate \
 --i-sequences fish-12S-ref-seqs_24Aug2023.qza \
 --i-taxa fish-12S-ref-tax_24Aug2023.qza \
 --p-mode 'uniq' \
 --p-threads 20 \
 --p-rank-handles 'disable' \
 --o-dereplicated-sequences fish-12S-ref-seqs_24Aug2023-derep.qza \
 --o-dereplicated-taxa fish-12S-ref-tax_24Aug2023-derep.qza
 
 #filter low-quality sequences and remove
 qiime rescript cull-seqs \
 --i-sequences fish-12S-ref-seqs_24Aug2023-derep.qza \
 --p-n-jobs 20 \
 --p-num-degenerates 5 \
 --p-homopolymer-length 8 \
 --o-clean-sequences fish-12S-ref-seqs_24Aug2023-cull.qza
 
 #now filter by sequence length
 qiime rescript filter-seqs-length \
 --i-sequences fish-12S-ref-seqs_24Aug2023-cull.qza \
 --p-global-min 90 \
 --p-global-max 1000 \
 --o-filtered-seqs fish-12S-ref-seqs_24Aug2023-FINAL.qza \
 --o-discarded-seqs fish-12S-ref-seqs_24Aug2023-discard.qza
 
#filter the derep taxonomy to include only the seqs in the FINAL ref file
qiime rescript filter-taxa \
--i-taxonomy fish-12S-ref-tax_24Aug2023-derep.qza \
--m-ids-to-keep-file fish-12S-ref-seqs_24Aug2023-FINAL.qza \
--o-filtered-taxonomy fish-12S-ref-tax_24Aug2023-FINAL.qza

#visualize for evaluation
qiime rescript evaluate-taxonomy \
--i-taxonomies fish-12S-ref-tax_24Aug2023-FINAL.qza \
--o-taxonomy-stats fish-12S-ref-tax_24Aug2023-FINAL-eval.qzv

#tabulate and visualize output
qiime metadata tabulate \
--m-input-file fish-12S-ref-tax_24Aug2023-FINAL.qza \
--o-visualization fish-12S-ref-tax_24Aug2023-FINAL.qzv &&
qiime rescript evaluate-seqs \
--i-sequences fish-12S-ref-seqs_24Aug2023-FINAL.qza \
--p-kmer-lengths 32 16 8 \
--o-visualization fish-12S-ref-seqs_24Aug2023-FINAL-eval.qzv
 


############################################
#############16S############################
############################################
qiime rescript get-ncbi-data \
    --p-query "txid7776[ORGN] AND (16S[Title] OR 16S ribosomal RNA[Title] OR 16S rRNA[Title]) AND (mitochondrion[Filter] OR plastid[Filter]) NOT environmental sample[Title] NOT environmental samples[Title] NOT environmental[Title] NOT uncultured[Title] NOT unclassified[Title] NOT unidentified[Title] NOT unverified[Title]" \
    --p-ranks kingdom phylum class order family genus species \
    --p-rank-propagation \
    --p-n-jobs 20 \
    --o-sequences fish-16S-ref-seqs.qza \
    --o-taxonomy fish-16S-ref-tax.qza \
    --verbose
    
#Dereplicate reference data
qiime rescript dereplicate \
 --i-sequences fish-16S-ref-seqs.qza \
 --i-taxa fish-16S-ref-tax.qza \
 --p-mode 'uniq' \
 --p-threads 20 \
 --p-rank-handles 'disable' \
 --o-dereplicated-sequences fish-16S-ref-seqs-derep.qza \
 --o-dereplicated-taxa fish-16S-ref-tax-derep.qza
 
 #filter low-quality sequences
 qiime rescript cull-seqs \
 --i-sequences fish-16S-ref-seqs-derep.qza \
 --p-n-jobs 20 \
 --p-num-degenerates 5 \
 --p-homopolymer-length 8 \
 --o-clean-sequences fish-16S-ref-seqs-cull.qza
 
 #now filter by sequence length
 qiime rescript filter-seqs-length \
 --i-sequences fish-16S-ref-seqs-cull.qza \
 --p-global-min 90 \
 --p-global-max 1000 \
 --o-filtered-seqs fish-16S-ref-seqs-keep.qza \
 --o-discarded-seqs fish-16S-ref-seqs-discard.qza
 
 qiime rescript evaluate-taxonomy \
 --i-taxonomies fish-16S-ref-seqs-keep.qza \
 --o-taxonomy-stats fish-16S-ref-seqs-keep-eval.qzv
 
 
############################################
#############COI############################
############################################
#Obtain a custom "fish" COI database using qiime and rescript
#txid77766 is Gnathostomata and includes sharks, rays, and all jawed fishes
#activate the rescript environment
conda activate rescript

qiime rescript get-ncbi-data \
    --p-query "txid7776[ORGN] AND\
    (COI[Title] OR COX[Title] OR COX1[Title] OR CO1[Title] OR\
    cytochrome c oxidase subunit 1[Title] OR cytochrome c oxidase subunit I[Title] OR\
    cytochrome oxidase subunit 1[Title] OR cytochrome oxidase subunit I[Title] OR\
    mitochondrion[Title] OR mitochondrial[Title]) AND\
    (mitochondrion[Filter] OR plastid[Filter] OR genomic dna[Filter]) NOT\
    environmental sample[Title] NOT\
    environmental samples[Title] NOT\
    environmental[Title] NOT\
    uncultured[Title] NOT\
    unclassified[Title] NOT\
    unidentified[Title] NOT\
    unverified[Title]"\
    --p-ranks kingdom phylum class order family genus species \
    --p-rank-propagation \
    --p-n-jobs 20 \
    --o-sequences fish-COI-ref-seqs-25Aug2023.qza \
    --o-taxonomy fish-COI-ref-tax-25Aug2023.qza \
    --verbose | at 22:00  #adding at command to run overnight

    
#Dereplicate reference data - sequences and taxonomy
qiime rescript dereplicate \
 --i-sequences fish-COI-ref-seqs-25Aug2023.qza \
 --i-taxa fish-COI-ref-tax-25Aug2023.qza \
 --p-mode 'uniq' \
 --p-threads 20 \
 --p-rank-handles 'disable' \
 --o-dereplicated-sequences fish-COI-ref-seqs-25Aug2023-derep.qza \
 --o-dereplicated-taxa fish-COI-ref-tax-25Aug2023-derep.qza
 
 #filter low-quality sequences and remove
 qiime rescript cull-seqs \
 --i-sequences fish-COI-ref-seqs-25Aug2023-derep.qza \
 --p-n-jobs 20 \
 --p-num-degenerates 5 \
 --p-homopolymer-length 8 \
 --o-clean-sequences fish-COI-ref-seqs-25Aug2023-cull.qza
 
 #now filter by sequence length
 qiime rescript filter-seqs-length \
 --i-sequences fish-COI-ref-seqs-25Aug2023-cull.qza \
 --p-global-min 90 \
 --p-global-max 1000 \
 --o-filtered-seqs fish-COI-ref-seqs-25Aug2023-FINAL.qza \
 --o-discarded-seqs fish-COI-ref-seqs-25Aug2023-discard.qza
 
#filter the derep taxonomy to include only the seqs in the FINAL ref file
qiime rescript filter-taxa \
--i-taxonomy fish-COI-ref-tax-25Aug2023-derep.qza \
--m-ids-to-keep-file fish-COI-ref-seqs-25Aug2023-FINAL.qza \
--p-exclude "unassigned" \
--o-filtered-taxonomy fish-COI-ref-tax-25Aug2023-FINAL.qza

#visualize for evaluation
qiime rescript evaluate-taxonomy \
--i-taxonomies fish-COI-ref-tax-25Aug2023-FINAL.qza \
--o-taxonomy-stats fish-COI-ref-tax-25Aug2023-FINAL-eval.qzv

qiime tools view fish-COI-ref-tax-25Aug2023-FINAL-eval.qzv

#evaluate
qiime rescript evaluate-seqs \
--i-sequences fish-COI-ref-seqs-25Aug2023-FINAL.qza \
--p-kmer-lengths 32 16 8 \
--o-visualization fish-COI-ref-seqs-25Aug2023-FINAL-eval.qzv

qiime tools view fish-COI-ref-seqs-25Aug2023-FINAL-eval.qzv


#tabulate and visualize output (in the regular qiime2 environment)
conda deactivate
conda activate qiime2_2023_5

qiime metadata tabulate \
--m-input-file fish-COI-ref-tax-25Aug2023-FINAL.qza \
--o-visualization fish-COI-ref-tax-25Aug2023-FINAL.qzv

qiime tools view fish-COI-ref-tax-25Aug2023-FINAL.qzv 

