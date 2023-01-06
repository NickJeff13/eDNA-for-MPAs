#removing primers from the 5' end of both forward and reverse reads as illumina adapters have already been removed
for i in *_R1.fastq.gz;
do
  SAMPLE=$(echo ${i} | sed "s/_R1.fastq.gz//") 
  echo ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz
  cutadapt -g ^AGCGYAATCACTTGTCTYTTAA -G ^CRBGGTCGCCCCAACCRAA \
    --discard-untrimmed --nextseq-trim=20 -m 20 -q 20,20 -n 2 \
    --match-read-wildcards \
     -o ${SAMPLE}_R1.fastq \
    -p ${SAMPLE}_R2.fastq \
    ${SAMPLE}_R1.fastq.gz  ${SAMPLE}_R2.fastq.gz ;
done

done
qiime cutadapt trim-paired \
--i-demultiplexed-sequences 16S-combined-demux.qza \
--p-cores 40 \
--p-anywhere-f AGCGYAATCACTTGTCTYTTAA \
--p-anywhere-r CRBGGTCGCCCCAACCRAA \
--p-error-rate 0.11 \
--p-discard-untrimmed True \
--p-match-read-wildcards \
--p-match-adapter-wildcards \
--o-trimmed-sequences 16s-demux-trimmed.qza \
--output-dir  trimmed \
--verbose
