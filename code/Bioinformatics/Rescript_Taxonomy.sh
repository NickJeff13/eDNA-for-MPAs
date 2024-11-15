############################################
#############COI############################
############################################
#Obtain a custom "fish" COI database using qiime and rescript
#txid77766 is Gnathostomata and includes sharks, rays, and all jawed fishes
#txid33208 is all metazoa
#txid50557 is insecta
#activate the rescript environment
conda activate rescript

qiime rescript get-ncbi-data \
    --p-query "txid33208[ORGN] NOT\
    txid50557[ORGN] AND\
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
    --p-n-jobs 40 \
    --o-sequences fish-COI-ref-seqs-19Sept2023.qza \
    --o-taxonomy fish-COI-ref-tax-19Sept023.qza \
    --verbose | at 22:00  #adding at command to run overnight

    
#Dereplicate reference data - sequences and taxonomy
qiime rescript dereplicate \
 --i-sequences fish-COI-ref-seqs-19Sept2023.qza \
 --i-taxa fish-COI-ref-tax-19Sept023.qza \
 --p-mode 'uniq' \
 --p-threads 40 \
 --p-rank-handles 'disable' \
 --o-dereplicated-sequences fish-COI-ref-seqs-19Sept2023-derep.qza \
 --o-dereplicated-taxa fish-COI-ref-tax-19Sept2023-derep.qza
 
#filter low-quality sequences and remove
 qiime rescript cull-seqs \
 --i-sequences fish-COI-ref-seqs-19Sept2023-derep.qza \
 --p-n-jobs 40 \
 --p-num-degenerates 5 \
 --p-homopolymer-length 8 \
 --o-clean-sequences fish-COI-ref-seqs-19Sept2023-culled.qza
 
#now filter by sequence length
 qiime rescript filter-seqs-length \
 --i-sequences fish-COI-ref-seqs-19Sept2023-culled.qza \
 --p-global-min 90 \
 --p-global-max 1000 \
 --o-filtered-seqs fish-COI-ref-seqs-19Sept2023-FINAL.qza \
 --o-discarded-seqs fish-COI-ref-seqs-19Sept2023-discard.qza
 
#filter the derep taxonomy to include only the seqs in the FINAL ref file
qiime rescript filter-taxa \
--i-taxonomy fish-COI-ref-tax-19Sept2023-derep.qza \
--m-ids-to-keep-file fish-COI-ref-seqs-19Sept2023-FINAL.qza \
--o-filtered-taxonomy fish-COI-ref-tax-19Sept2023-FINAL.qza

#export the fasta and taxonomy tsv
qiime tools export \
--input-path fish-COI-ref-seqs-19Sept2023-FINAL.qza \
--output-path fish-COI-ref-seqs-19Sept2023-FINAL.fasta

qiime tools export \
--input-path fish-COI-ref-tax-19Sept2023-FINAL.qza \
--output-path fish-COI-ref-tax-19Sept2023-FINAL.tsv

#remove line breaks in the middle of sequences using awk
awk '/^[>;]/ { if (seq) { print seq }; seq=""; print } /^[^>;]/ { seq = seq $0 } END { print seq }' fish-COI-ref-seqs-19Sept2023-FINAL.fasta > fish-COI-ref-seqs-19Sept2023-FINAL-clean.fasta

#visualize for evaluation
qiime rescript evaluate-taxonomy \
--i-taxonomies fish-COI-ref-tax-19Sept2023-FINAL.qza \
--o-taxonomy-stats fish-COI-ref-tax-19Sept2023-FINAL-eval.qzv

qiime tools view fish-COI-ref-tax-19Sept2023-FINAL-eval.qzv

#evaluate
qiime rescript evaluate-seqs \
--i-sequences fish-COI-ref-seqs-19Sept2023-FINAL.qza \
--p-kmer-lengths 32 16 8 \
--o-visualization fish-COI-ref-seqs-19Sept2023-FINAL-eval.qzv

qiime tools view fish-COI-ref-seqs-19Sept2023-FINAL-eval.qzv

#tabulate and visualize output
qiime metadata tabulate \
--m-input-file fish-COI-ref-tax-19Sept2023-FINAL.qza \
--o-visualization fish-COI-ref-tax-19Sept2023-FINAL.qzv

qiime tools view fish-COI-ref-tax-19Sept2023-FINAL.qzv 

