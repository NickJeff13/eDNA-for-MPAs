#!/bin/bash
#Use cutadapt to trim primers for 12S MiFish
#MiFish forward NNNNNNGTCGGTAAAACTCGTGCCAGC
#MiFish reverse NNNNNNCATAGTGGGGTATCTAATCCCAGTTTG
:'
If you have multiple underscores in your file names, but the number of underscores is consistent, you could simply tweak the number given to cut -f. 
Otherwise, you might need to use sed to remove parts of the names, or even use rename to batch rename the files to something that will work
'
rm -r output
samples=$(ls | cut -dR -f1 | sort | uniq)
echo $samples

#run cutadapt in loop - note the '^' in the primer sequence means it was anchored to the adapter(s) 
mkdir output
for s in $samples;
do
        cutadapt -g ^NNNNNNGTCGGTAAAACTCGTGCCAGC -G ^NNNNNNCATAGTGGGGTATCTAATCCCAGTTTG \
        -o output/${s}R1.fastq -p output/${s}R2.fastq --discard-untrimmed \
        ${s}R1.fastq ${s}R2.fastq ;
done


:'
Example output:
=== First read: Adapter 1 ===

Sequence: NNNNNNGTCGGTAAAACTCGTGCCAGC; Type: anchored 5; Length: 27; Trimmed: 535031 times

No. of allowed errors: 2

Overview of removed sequences
length	count	expect	max.err	error counts
25	239	0.0	2	0 0 239
26	5027	0.0	2	0 4665 362
27	528486	0.0	2	499820 26529 2137
28	1267	0.0	2	0 1127 140
29	12	0.0	2	0 0 12


=== Second read: Adapter 2 ===

Sequence: NNNNNNCATAGTGGGGTATCTAATCCCAGTTTG; Type: anchored 5; Length: 33; Trimmed: 525079 times

No. of allowed errors: 2

Overview of removed sequences
length	count	expect	max.err	error counts
31	475	0.0	2	0 0 475
32	7550	0.0	2	0 7315 235
33	515191	0.0	2	499655 14786 750
34	1855	0.0	2	0 1781 74
35	8	0.0	2	0 0 8
' 
#!/usr/bin/bash
#removing primers from the 5' end of both forward and reverse reads as illumina adapters have already been removed
#Can tailor this to remove adapters or just one primer if only using forward reads
for i in *_R1.fastq.gz;
do
  SAMPLE=$(echo ${i} | sed "s/_R1.fastq.gz//") 
  echo ${SAMPLE}_R1.fastq.gz ${SAMPLE}_R2.fastq.gz
  cutadapt -g GTCGGTAAAACTCGTGCCAGC \
  -G NNNNNNCATAGTGGGGTATCTAATCCCAGTTTG \
    --discard-untrimmed --nextseq-trim 20 -m 20 -n 2 \
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
